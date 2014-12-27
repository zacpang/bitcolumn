//
//  print_helper.cpp
//  bitcolumn
//
//  Created by ZhifeiPang on 12/15/14.
//  Copyright (c) 2014 ZhifeiPang. All rights reserved.
//

#include "print_helper.h"


/**
 *  Pirnt 1-d array
 **/
void print(uint32_t *a, int count)
{
    for (int i = 0; i<count; i++)
    {
        if (i % 32 == 0)
            printf("\n------------------------------------\n");
        printf("%u\t", a[i]);
    }
}

/**
 *  Print 2-d array
 **/
void print2d(uint32_t a[][N_BITS], int row_count)
{
    for (int i=0; i<row_count; i++)
    {
        for (int j = 0; j<N_BITS; j++)
        {
            printf("%u\t", a[i][j]);
        }
        printf("\n");
    }
}


void print_array_data(uint8_t array_data[N][CACHE_LINE_SIZE])
{
	printf("array_data:\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N_ATTR; j++)
		{
			for (int k = 0; k < N_BYTE_INARRAY; k++)
				printf("%u ", array_data[i][j*N_BYTE_INARRAY + k]);
			printf("  ");
		}
		printf("\n");
	}
}

/**
 * Print 2-d simd matrix
 **/
void print3dmatrix(__m128i matrix[N_ATTR][N_ROWS][N_BITS])
{
	for (int m = 0; m < N_ATTR; m++)
	{
		printf("\nAttr:%d\n",m);
		for (int k = 0; k < N_ROWS; k++)
		{
			for (int i = 0; i < N_SLOTS; i++)
			{
				for (int j = 0; j < N_BITS; j++)
				{
					uint32_t *t = (uint32_t*)&matrix[m][k][j];
					printf("%u,", t[i]);
				}
				printf("\n\n\n");
			}
			printf("********************************ROW**END*********************************\n");
		}
	}
}

void print2dmatrix(__m128i v[N_ATTR][N_BITS])
{
	for (int i = 0; i < N_ATTR; i++)
	{
		printf("attr2d%d\n", i);
		for (int j = 0; j < N_BITS; j++)
		{
			uint32_t *t = (uint32_t*)&v[i][j];
			for (int k = 0; k < N_SLOTS; k++)
			{
				printf("%u\t", t[i]);
			}
		}
	}
}

/**
 * Print 1-d simd matrix
 **/
void print1dmatrix(__m128i matrix[N_ROWS])
{
	for (int i = 0; i<N_ROWS; i++)
    {
        for (int j=0; j<N_SLOTS; j++)
        {
            uint32_t *t = (uint32_t*)&matrix[i];
            printf("%u\t", t[j]);
        }
    }
    printf("\n");
}
