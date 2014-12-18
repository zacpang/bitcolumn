//
//  Header.h
//  bitcolumn
//
//  Created by ZhifeiPang on 12/15/14.
//  Copyright (c) 2014 ZhifeiPang. All rights reserved.
//

#ifndef bitcolumn_Header_h
#define bitcolumn_Header_h

#include "bitcolumn.h"
#include <xmmintrin.h>
#include <inttypes.h>
#include <stdio.h>


void print(uint32_t *a, int count);
void print2d(uint32_t a[][N_BITS], int row_count);
void print2dmatrix(__m128i matrix[][N_BITS]);
void print1dmatrix(__m128i matrix[], int len);


#endif
