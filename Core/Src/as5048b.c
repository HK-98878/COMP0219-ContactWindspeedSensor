#include "stm32f4xx_hal.h"
#include "stm32f4xx_hal_def.h"
#include <stdint.h>

#define AS5048B_AGC_ADDR 0xFA
#define AS5048B_DIAG_ADDR 0xFB
#define AS5048B_ANGLE_ADDR 0xFE

#define AS5048B_ANGLE_DATA_SIZE 2
#define AS5048B_SAFE_ANGLE_DATA_SIZE (AS5048B_ANGLE_ADDR - AS5048B_DIAG_ADDR + AS5048B_ANGLE_DATA_SIZE)
#define AS5048B_AGC_ANGLE_DATA_SIZE (AS5048B_ANGLE_ADDR - AS5048B_AGC_ADDR + AS5048B_ANGLE_DATA_SIZE)

#define EXTRACT_ANGLE_DATA(data, idx) (((data)[(idx)] << 6) | ((data)[(idx)+1] & 0x3F))


HAL_StatusTypeDef AS5048_CheckSetup(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr) {
  uint8_t data;
  HAL_StatusTypeDef status = HAL_I2C_Mem_Read(i2c_handle, sensor_addr, AS5048B_DIAG_ADDR,
    I2C_MEMADD_SIZE_8BIT, &data, 1, 100);
  if(status == HAL_OK && data & 0b1) return HAL_OK;
  return HAL_BUSY;
}

uint16_t AS5048_ConvertAngleBuffer(uint8_t angle_buffer[2]) {
  return EXTRACT_ANGLE_DATA(angle_buffer, 0);
}
HAL_StatusTypeDef AS5048_ReadAngle(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint16_t *angle) {
  uint8_t data[AS5048B_ANGLE_DATA_SIZE]; 
  HAL_StatusTypeDef status = HAL_I2C_Mem_Read(i2c_handle, sensor_addr, AS5048B_ANGLE_ADDR,
    I2C_MEMADD_SIZE_8BIT, data, AS5048B_ANGLE_DATA_SIZE, 100);
  if (status == HAL_OK){
    *angle = EXTRACT_ANGLE_DATA(data, 0);
  }
  return status;
}
HAL_StatusTypeDef AS5048_ReadAngleDMA(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint8_t angle_buffer[2]) {
  return HAL_I2C_Mem_Read_DMA(i2c_handle, sensor_addr, AS5048B_ANGLE_ADDR, 
    I2C_MEMADD_SIZE_8BIT, angle_buffer, AS5048B_ANGLE_DATA_SIZE);
}

HAL_StatusTypeDef AS5048_ReadAngle_Safe(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint16_t *angle) {
  uint8_t data[AS5048B_SAFE_ANGLE_DATA_SIZE]; 
  HAL_StatusTypeDef status = HAL_I2C_Mem_Read(
    i2c_handle, sensor_addr, AS5048B_DIAG_ADDR,
    I2C_MEMADD_SIZE_8BIT, data, AS5048B_SAFE_ANGLE_DATA_SIZE, 100);
  if (status == HAL_OK){
    if(data[0] & 0b11) {
      *angle = 0xffff;
      return HAL_ERROR;
    }
    *angle = EXTRACT_ANGLE_DATA(data, AS5048B_ANGLE_ADDR-AS5048B_DIAG_ADDR);
  }
  return status;
}

HAL_StatusTypeDef AS5048_ReadAngle_AGC(I2C_HandleTypeDef* i2c_handle, const uint8_t sensor_addr, uint16_t *angle, uint8_t *agc) {
  uint8_t data[AS5048B_AGC_ANGLE_DATA_SIZE]; 
  HAL_StatusTypeDef status = HAL_I2C_Mem_Read(
    i2c_handle, sensor_addr, AS5048B_AGC_ADDR,
    I2C_MEMADD_SIZE_8BIT, data, AS5048B_AGC_ANGLE_DATA_SIZE, 100);
  if (status == HAL_OK){
    *agc = data[0];
    *angle = EXTRACT_ANGLE_DATA(data, AS5048B_ANGLE_ADDR-AS5048B_AGC_ADDR);
    // (data[AS5048B_ANGLE_ADDR-AS5048B_AGC_ADDR] << 6) | (data[AS5048B_ANGLE_ADDR-AS5048B_AGC_ADDR+1] & 0x3F);
  }
  return status;
}