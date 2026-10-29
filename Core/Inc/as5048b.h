#ifndef __as5048b_driver_h
#define __as5048b_driver_h

#include "stm32f4xx_hal.h"
#include "stm32f4xx_hal_def.h"
#include <stdint.h>

HAL_StatusTypeDef AS5048_CheckSetup(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr);
uint16_t AS5048_ConvertAngleBuffer(uint8_t angle_buffer[2]);

/* Read the 14-bit angle from the sensor
 *
 * Fastest possible version (compared to AS5048_ReadAngle_Safe and AS5048_ReadAngle_AGC), at the cost of safety checks
 *
 * Reads a 2 byte block to ensure matched reading bytes (no interleaved sensor update)
 */
HAL_StatusTypeDef AS5048_ReadAngle(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint16_t *angle);
HAL_StatusTypeDef AS5048_ReadAngleDMA(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint8_t angle_buffer[2]);

/* Read the 14-bit angle from the sensor, with a status check
 *
 * Returns HAL_ERROR if either the Comp High or Comp Low flags are set, indicating a weak or over-strong magnetic 
 * field respectively. Angle value is not updated in this case, as it cannot be trusted. To differentiate from an 
 * I2C-raised error, angle is also set to 0xffff (out of range for possible sensor outputs)
 *
 * Reads a 5 byte block to ensure matched readings+status flags (no interleaved sensor update)
 * - Also faster than segmented reads for reading and status due to addressing overhead
 */
HAL_StatusTypeDef AS5048_ReadAngle_Safe(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint16_t *angle);

/* Read the 14-bit angle and the AGC value from the sensor
 *
 * AGC (automatic gain control) provides an indication of magnetic field strength to inform reading reliability
 *
 * Reads a 6-byte block to ensure matched readings+status flags (no interleaved sensor update)
 * - Also technically fractionally faster than segmented read - same net clock pulses 
 *   (given the data spacing in memory), plus fewer START/STOP bit signals
 */
HAL_StatusTypeDef AS5048_ReadAngle_AGC(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint16_t *angle, uint8_t *agc);


#endif