//
//  bitcolumn.h
//  bitcolumn
//
//  Created by ZhifeiPang on 12/15/14.
//  Copyright (c) 2014 ZhifeiPang. All rights reserved.
//

#ifndef bitcolumn_bitcolumn_h
#define bitcolumn_bitcolumn_h

#include <stdio.h>
#include <inttypes.h>
#include <memory.h>
#include <xmmintrin.h>
#include <time.h>
#include <stdlib.h>
#include <emmintrin.h>


#define _2POW31_ 2147483648

#define N_BITS 32                                  //number of bits of uint32_t
#define N_SLOTS (128/NB)
#define UINT32MAX 4294967295                         //the MAX value of uint32_t
#define V_LEN 128
#define N_ROWS ((N + 127) / 128)                   //the number of raws of simd matrix
#define N 128                                      //the # of numbers
#define N_ATTR  1                                  // The number of the attributes


#endif
