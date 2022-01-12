/*******************************************************************************
 * Copyright (C) 2020 Mariusz Matusiak
 * @file    buffers.h - header file
 * @author  Mariusz Matusiak <mmatusiak@iis.p.lodz.pl>
 *                                     ORCID: 0000-0003-2407-3135
 * @brief   buffers.h header file consisting of the declarations of functions
 *          for memory management, allocating new and reallocating the used
 *          resources
 * @date    2018/09/24
 ******************************************************************************/
#ifndef INC_BUFFERS_H_
#define INC_BUFFERS_H_

#include "types.h"

#ifdef STM32F7
#define TOTAL_MCU_SRAM_IN_KB 320
#elif defined STM32H7
#define TOTAL_MCU_SRAM_IN_KB 864
#else
#define TOTAL_MCU_SRAM_IN_KB 320
#endif
#define TOTAL_MCU_SRAM    (TOTAL_MCU_SRAM_IN_KB << 10)
#define SRAM_ALLOC_STACK  20u
#define WARNING_RAM_DIV   1  /* TOTAL_MCU_SRAM >> (x) (1)-/2 (2)-/4 (3)-/8 */
extern const size_t compile_SRAM_allocated;  //42912

struct memory_record {
    uint32_t ptr_address;
    size_t   bytes;
};
typedef struct memory_record mem_record;


void * memory_allocator  (const size_t number_of_elements, const size_t
        size_of_element_in_bytes);
void   memory_deallocator(void * ptr);

#endif /* INC_BUFFERS_H_ */
