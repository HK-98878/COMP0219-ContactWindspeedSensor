/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2026 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "as5048b.h"
#include "stm32f4xx_hal_def.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
I2C_HandleTypeDef hi2c1;
DMA_HandleTypeDef hdma_i2c1_rx;

TIM_HandleTypeDef htim3;

UART_HandleTypeDef huart2;

/* USER CODE BEGIN PV */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_I2C1_Init(void);
static void MX_USART2_UART_Init(void);
static void MX_TIM3_Init(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

void UART_Print(const char *message) {
  HAL_UART_Transmit(&huart2 , (uint8_t *)message,
    (uint16_t)strlen(message), 1000);
}
#define ENCODER_I2C_ADDR (0b1000000 << 1) // 1 bit Lshift for HAL API
//
// I2C ISRs for non-blocking DMA read
//
static uint8_t rx_buf[2];
static uint32_t rx_buf_seq_idx = 0;

static volatile uint8_t  i2c_busy = 0;
static volatile uint32_t overruns = 0, i2c_errors = 0;

// Block accumulation state. Shared between the TIM3 and I2C1 ISRs
// Skipped/failed reads still "accumulate" in timer ticks, so block is always dispatched at
// 200Hz, even if I2C delays/fails
static struct accumulation_block cur;
static uint8_t  tick_in_block = 0;  // 0..AVG_N-1, position of the next tick in the block
static uint32_t block_id = 0;       // incremented every time a block is closed
static uint32_t inflight_block_id;  // block_id at the time the in-flight read was started

static void CloseBlock(void) {
  if (cur.n > 0) {                  // all reads failed -> no output for this period
    if (!RingPush(&cur)) ++overruns;
  }
  cur.n = 0;
  cur.sum = 0;
  cur.start_idx = rx_buf_seq_idx;
  ++block_id;
}

void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim) {
  if (htim->Instance != TIM3) return;
  rx_buf_seq_idx += 1; // Increment this regardless of success/fail on the read to detect skipped reads

  // Every AVG_N ticks, the previous block has had all its reads issued (and, at
  // ~120 us per read vs a 500 us period, completed), so close it now.
  if (tick_in_block == 0) CloseBlock();
  if (++tick_in_block == AVG_N) tick_in_block = 0;

  if (i2c_busy) { overruns++; return; } // previous read not finished
  i2c_busy = 1;
  inflight_block_id = block_id;
  if (AS5048_ReadAngleDMA(&hi2c1, ENCODER_I2C_ADDR, rx_buf) != HAL_OK) {
    i2c_busy = 0;
    overruns++;
  }
}

void HAL_I2C_MemRxCpltCallback(I2C_HandleTypeDef *hi2c) {
  static uint16_t enc_ref;
  static uint8_t enc_ref_set = 0;
  if (hi2c->Instance != I2C1) return;
  i2c_busy = 0;

  // A read that was still in flight when its block closed belongs to a block
  // that has already been pushed; drop it rather than mis-attribute it.
  if (inflight_block_id != block_id) { overruns++; return; }

  uint16_t raw = AS5048_ConvertAngleBuffer(rx_buf);

  // Fixed reference, latched once from the first sample after boot, used 
  // for wrap detection/unwrapping
  if(!enc_ref_set) {
    enc_ref = raw;
    enc_ref_set = 1;
  }
  cur.ref = enc_ref;

  // Manage encoder wrap (modular difference to the fixed ref, signed result)
  int32_t d = (int32_t)raw - cur.ref;
  if(d > ENC_RANGE / 2) d -= ENC_RANGE;
  if(d < -ENC_RANGE / 2) d += ENC_RANGE;
  cur.sum += d;
  ++cur.n;
}

void HAL_I2C_ErrorCallback(I2C_HandleTypeDef *hi2c) {
  if (hi2c->Instance != I2C1) return;
  i2c_errors++;
  i2c_busy = 0;   // consider bus recovery here if errors persist
}

// Air density calibrated in the VDAS
#define RHO 1.2056f
#ifndef RHO // Fallback
#define RHO 1.2056f      // kg/m^3, VDAS reference density
#endif

// Constants
#define SQRT_A 0.059487f  // Pa/count

#define QUAD_A1  0.0447421f
#define QUAD_A2  2.12011e-5f

// v = k sqrt(theta), theta = c0 - filtered count (counts, positive with airspeed)
float sqrt_model(float theta) {
  const static float c = 2.0f * SQRT_A / RHO;

  if(theta <= 0.0f) return 0.0f;
  return sqrtf(c * theta);
}

// v = k sqrt(theta) sqrt(1 + b theta), theta = c0 - filtered count
float quadratic_model(float theta) {
  // Roll k inside the sqrt for fewer calculations
  const static float c1 = 2.0f * QUAD_A1 / RHO;
  const static float c2 = 2.0f * QUAD_A2 / RHO;

  if(theta <= 0.0f) return 0.0f;
  return sqrtf(theta * (c1 + c2 * theta));
}

float cubic_model(float theta) {
  const static float scale = 282.1641f;
  const static float a = 1.37421e-8 / scale;
  const static float b = -0.00043698 / scale;
  const static float c = 4.63856 / scale;

  if(theta <= 0.0f) return 0.0f;
  return a*powf(theta, 3.0) + b*powf(theta, 2.0) + c*theta;
}



//
// Tare (c0) procedure
//
// c0 = mean count over a still-air period, taken at boot and on each B1 press.
// Never re-zeroed automatically. Output is silent while taring.
//
#define TARE_BLOCKS        400   // tare period in blocks: 400 x AVG_N samples ~= 2 s
#define TARE_BLINK_BLOCKS  20    // LD2 toggle period while taring (~100 ms)
// Stillness check: std dev of the 10-sample block means over the tare period.
#define TARE_MAX_STD       2.0f  // counts

struct tare_state {
  int32_t  sum;      // sum of block sums (counts rel. to fixed ref); max 400*10*8192 ~= 33M
  uint32_t n;        // samples accumulated
  uint16_t blocks;   // blocks accumulated
  float    mean;     // Welford running mean of block means (stillness check only)
  float    m2;       // Welford sum of squared deviations
};

static void TareReset(struct tare_state *t) {
  memset(t, 0, sizeof(*t));
}

// Accumulate one block. Returns 1 when the tare period is complete and still
// (c0 written to *c0_off), 0 otherwise. A failed period resets and retries.
static int TareAccumulate(struct tare_state *t, const struct accumulation_block *b, float *c0_off) {
  t->sum += b->sum;
  t->n   += b->n;
  ++t->blocks;

  float m = (float)b->sum / b->n;
  float delta = m - t->mean;
  t->mean += delta / t->blocks;
  t->m2   += delta * (m - t->mean);

  if (t->blocks % TARE_BLINK_BLOCKS == 0) HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
  if (t->blocks < TARE_BLOCKS) return 0;

  float var = t->m2 / (t->blocks - 1);
  int still = var <= TARE_MAX_STD * TARE_MAX_STD;
  if (still) {
    // double for the division: sum can exceed float's 24-bit exact range
    *c0_off = (float)((double)t->sum / (double)t->n);
  }
  TareReset(t);
  return still;
}

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */
  HAL_DBGMCU_EnableDBGSleepMode(); // Prevents board from locking under debug mode due to __WFI
  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_I2C1_Init();
  MX_USART2_UART_Init();
  MX_TIM3_Init();
  /* USER CODE BEGIN 2 */

  // Check encoder is connected and active. No UART messages: the downstream
  // device only accepts single float lines. LD2 stays off while stuck here.
  HAL_StatusTypeDef status;
  status = HAL_I2C_IsDeviceReady (&hi2c1, ENCODER_I2C_ADDR, 3, 100);
  if (status != HAL_OK) {
    while (1);
  }
  do {
    status = AS5048_CheckSetup(&hi2c1, ENCODER_I2C_ADDR);
  } while(status != HAL_OK);

  char msg[100];

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  
  struct accumulation_block data;
  struct moving_window window = {0};

  // Tare runs first (boot), and again on each B1 press. The tare period (400
  // blocks) is far longer than the moving window (10 blocks), so the window is
  // always full by the time output starts - no separate warm-up state needed.
  enum { STATE_TARING, STATE_RUNNING } state = STATE_TARING;
  struct tare_state tare;
  TareReset(&tare);
  float c0_off = 0.0f;              // c0 as an offset from the fixed encoder ref
  uint8_t b1_prev = 0;

  HAL_TIM_Base_Start_IT(&htim3); // Start sensor read timer

  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
    while(RingPop(&data)) {
      MovingWindowAccumulate(&window, &data);

      // B1 (active low) polled once per block (5 ms). Press -> restart tare
      uint8_t b1 = HAL_GPIO_ReadPin(B1_GPIO_Port, B1_Pin) == GPIO_PIN_RESET;
      if (b1 && !b1_prev) {
        TareReset(&tare);
        HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);
        state = STATE_TARING;
      }
      b1_prev = b1;

      if (state == STATE_TARING) {
        if (TareAccumulate(&tare, &data, &c0_off)) {
          HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_SET);
          state = STATE_RUNNING;
        }
        continue;                   // silent while taring
      }

      // Filtered count as an offset from the fixed ref (exact integer sums,
      // single float conversion). n_total > 0: CloseBlock never pushes an
      // empty block. theta = c0 - count, positive with airspeed.
      float filtered_off = (float)window.sum_total / window.n_total;
      float v = quadratic_model(c0_off - filtered_off);
      // float v = sqrt_model(c0_off - filtered_off);
      // float v = cubic_model(c0_off - filtered_off);

      // This MUST be in main loop for heap allocation - newlib not reentrant
      snprintf(msg, sizeof(msg), "%.2f\r\n", (double)v);
      UART_Print(msg);
    }

    __WFI(); // Wait for any interrupt (after ISR completion)
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_BYPASS;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 8;
  RCC_OscInitStruct.PLL.PLLN = 336;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV4;
  RCC_OscInitStruct.PLL.PLLQ = 4;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief I2C1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_I2C1_Init(void)
{

  /* USER CODE BEGIN I2C1_Init 0 */

  /* USER CODE END I2C1_Init 0 */

  /* USER CODE BEGIN I2C1_Init 1 */

  /* USER CODE END I2C1_Init 1 */
  hi2c1.Instance = I2C1;
  hi2c1.Init.ClockSpeed = 400000;
  hi2c1.Init.DutyCycle = I2C_DUTYCYCLE_2;
  hi2c1.Init.OwnAddress1 = 0;
  hi2c1.Init.AddressingMode = I2C_ADDRESSINGMODE_7BIT;
  hi2c1.Init.DualAddressMode = I2C_DUALADDRESS_DISABLE;
  hi2c1.Init.OwnAddress2 = 0;
  hi2c1.Init.GeneralCallMode = I2C_GENERALCALL_DISABLE;
  hi2c1.Init.NoStretchMode = I2C_NOSTRETCH_DISABLE;
  if (HAL_I2C_Init(&hi2c1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN I2C1_Init 2 */

  /* USER CODE END I2C1_Init 2 */

}

/**
  * @brief TIM3 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM3_Init(void)
{

  /* USER CODE BEGIN TIM3_Init 0 */

  /* USER CODE END TIM3_Init 0 */

  TIM_ClockConfigTypeDef sClockSourceConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM3_Init 1 */

  /* USER CODE END TIM3_Init 1 */
  htim3.Instance = TIM3;
  htim3.Init.Prescaler = 41;
  htim3.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim3.Init.Period = 999;
  htim3.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim3.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim3) != HAL_OK)
  {
    Error_Handler();
  }
  sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
  if (HAL_TIM_ConfigClockSource(&htim3, &sClockSourceConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim3, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM3_Init 2 */

  /* USER CODE END TIM3_Init 2 */

}

/**
  * @brief USART2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART2_UART_Init(void)
{

  /* USER CODE BEGIN USART2_Init 0 */

  /* USER CODE END USART2_Init 0 */

  /* USER CODE BEGIN USART2_Init 1 */

  /* USER CODE END USART2_Init 1 */
  huart2.Instance = USART2;
  huart2.Init.BaudRate = 115200;
  huart2.Init.WordLength = UART_WORDLENGTH_8B;
  huart2.Init.StopBits = UART_STOPBITS_1;
  huart2.Init.Parity = UART_PARITY_NONE;
  huart2.Init.Mode = UART_MODE_TX_RX;
  huart2.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart2.Init.OverSampling = UART_OVERSAMPLING_16;
  if (HAL_UART_Init(&huart2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART2_Init 2 */

  /* USER CODE END USART2_Init 2 */

}

/**
  * Enable DMA controller clock
  */
static void MX_DMA_Init(void)
{

  /* DMA controller clock enable */
  __HAL_RCC_DMA1_CLK_ENABLE();

  /* DMA interrupt init */
  /* DMA1_Stream0_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream0_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream0_IRQn);

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};
  /* USER CODE BEGIN MX_GPIO_Init_1 */

  /* USER CODE END MX_GPIO_Init_1 */

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : B1_Pin */
  GPIO_InitStruct.Pin = B1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(B1_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LD2_Pin */
  GPIO_InitStruct.Pin = LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LD2_GPIO_Port, &GPIO_InitStruct);

  /* USER CODE BEGIN MX_GPIO_Init_2 */

  /* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
