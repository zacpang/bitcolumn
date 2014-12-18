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



#define N 256                                      //the # of data

#define _2POW31_ 2147483648
#define DATATYPE_LEN 32 
#define N_BITS 32                                  //含义已经变化。表示每一个ＳＩＭＤ矩阵的列数
#define N_SLOTS 4                                  //N_SLOTS*N_BITS MUST BE 128
#define UINT32MAX 4294967295                       //the MAX value of uint32_t
#define V_LEN 128
#define N_ROWS ((N + 127) / 128)                   //the number of raws of simd matrix

#define N_ATTR  1                                  // The number of the attributes

#endif
