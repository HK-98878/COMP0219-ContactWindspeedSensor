#include "main.h"

// 
// Ring buffer used to stash values from the I2C ISR
//
#define RING_BUFFER_COUNT 64
static struct {
  struct accumulation_block buffer[RING_BUFFER_COUNT]; // Buffering raw encoder counts
  volatile uint16_t head;
  volatile uint16_t tail;
} ring_buffer;
int RingPush(struct accumulation_block *data) {
  uint16_t next = (ring_buffer.head + 1) % RING_BUFFER_COUNT;
  if (next == ring_buffer.tail) return 0;      // full
  ring_buffer.buffer[ring_buffer.head] = *data;
  __DMB();
  ring_buffer.head = next;
  return 1;
}
int RingPop(struct accumulation_block* data) {
  if (ring_buffer.tail == ring_buffer.head) return 0;  // empty
  *data = ring_buffer.buffer[ring_buffer.tail];
  __DMB();
  ring_buffer.tail = (ring_buffer.tail + 1) % RING_BUFFER_COUNT;
  return 1;
}

//
// Moving window accumulation function for faster averaging
//

void MovingWindowAccumulate(struct moving_window* window, const struct accumulation_block* block) {
  if(window->blocks == MOVING_WINDOW_N) {
    // Evict the oldest block
    window->sum_total -= window->block_sum[window->idx];
    window->n_total   -= window->block_n[window->idx];
  } else {
    ++window->blocks;
  }

  window->block_sum[window->idx] = block->sum;
  window->block_n[window->idx]   = block->n;
  window->sum_total += block->sum;
  window->n_total   += block->n;
  window->idx = (window->idx + 1) % MOVING_WINDOW_N;
}