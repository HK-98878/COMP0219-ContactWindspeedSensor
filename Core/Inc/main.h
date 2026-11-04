/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.h
  * @brief          : Header for main.c file.
  *                   This file contains the common defines of the application.
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

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __MAIN_H
#define __MAIN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Includes ------------------------------------------------------------------*/
#include "stm32f4xx_hal.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Exported types ------------------------------------------------------------*/
/* USER CODE BEGIN ET */
#define MOVING_WINDOW_N 10
#define AVG_N 10
#define ENC_RANGE 16384

struct accumulation_block {
  uint32_t start_idx;
  int32_t sum;
  uint16_t ref;
  uint8_t n;
};

// Trailing moving window over the last MOVING_WINDOW_N blocks. Holds the raw
// integer block sums and sample counts, so the running totals are exact (no
// float drift) and every sample is weighted equally even when a block is short.
// Max |sum_total| = MOVING_WINDOW_N * AVG_N * ENC_RANGE/2 = 819200 -> int32.
struct moving_window {
  int32_t block_sum[MOVING_WINDOW_N];
  uint8_t block_n[MOVING_WINDOW_N];
  uint8_t blocks;      // number of blocks currently in the window (fills to MOVING_WINDOW_N)
  uint8_t idx;         // slot the next block will be written to
  int32_t sum_total;   // sum of block_sum[] (counts relative to the block ref)
  uint16_t n_total;    // sum of block_n[] (samples in the window)
};
/* USER CODE END ET */

/* Exported constants --------------------------------------------------------*/
/* USER CODE BEGIN EC */

/* USER CODE END EC */

/* Exported macro ------------------------------------------------------------*/
/* USER CODE BEGIN EM */

/* USER CODE END EM */

/* Exported functions prototypes ---------------------------------------------*/
void Error_Handler(void);

/* USER CODE BEGIN EFP */

int RingPush(struct accumulation_block *data);
int RingPop(struct accumulation_block* data);
void MovingWindowAccumulate(struct moving_window* window, const struct accumulation_block* block);

/* USER CODE END EFP */

/* Private defines -----------------------------------------------------------*/
#define B1_Pin GPIO_PIN_13
#define B1_GPIO_Port GPIOC
#define USART_TX_Pin GPIO_PIN_2
#define USART_TX_GPIO_Port GPIOA
#define USART_RX_Pin GPIO_PIN_3
#define USART_RX_GPIO_Port GPIOA
#define LD2_Pin GPIO_PIN_5
#define LD2_GPIO_Port GPIOA
#define TMS_Pin GPIO_PIN_13
#define TMS_GPIO_Port GPIOA
#define TCK_Pin GPIO_PIN_14
#define TCK_GPIO_Port GPIOA
#define SWO_Pin GPIO_PIN_3
#define SWO_GPIO_Port GPIOB

/* USER CODE BEGIN Private defines */

/* USER CODE END Private defines */

#ifdef __cplusplus
}
#endif

#endif /* __MAIN_H */
