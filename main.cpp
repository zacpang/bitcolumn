#include <stdio.h>
#include <inttypes.h>
#include <memory.h>
#include <xmmintrin.h>
#include <time.h>
#include <stdlib.h>

#define NB 32                                  //number of bits of uint32_t
#define SLOT (128/NB)
#define MAX 4294967295                         //the MAX value of uint32_t
#define N 51200000                           //the # of number
#define NR ((N + 127) / 128)                   //the number of raws of simd matrix

static uint32_t mask[NB] = {
	0x80000000, 0x40000000, 0x20000000, 0x10000000,
	0x08000000, 0x04000000, 0x02000000, 0x01000000,
	0x00800000, 0x00400000, 0x00200000, 0x00100000,
	0x00080000, 0x00040000, 0x00020000, 0x00010000,
	0x00008000, 0x00004000, 0x00002000, 0x00001000,
	0x00000800, 0x00000400, 0x00000200, 0x00000100,
	0x00000080, 0x00000040, 0x00000020, 0x00000010,
	0x00000008, 0x00000004, 0x00000002, 0x00000001,
};




__m128i res[NR];								//保存运用SSE时，比较的结果
uint32_t res_without_sse[N];                    //保存顺序查询时，比较的结果
__m128i matrix[NR][NB]; 
__m128i com_value[NB];                         //由value得到的，用于比较的数组
uint32_t data[N];
FILE *readdata;
/**
*  Pirnt 1-d array
**/
void print(uint32_t *a, int count)
{
	int i;
	for (i = 0; i<count; i++)
	{
		if (i % 128 == 0)
			printf("\n------------------------------------\n");
		printf("%u\t", a[i]);
	}
	printf("\n\n\n");

}

/**
*  Print 2-d array
**/
void print2d(uint32_t a[][32], int d)
{
	int i, j;
	for (i = 0; i<d; i++)
	{
		for (j = 0; j<32; j++)
		{
			printf("%u\t", a[i][j]);
		}
		printf("\n");
	}
	printf("\n\n\n");
}


/**
* Print 2-d simd matrix
**/
void print2dmatrix(__m128i matrix[][NB])
{

	for (int k = 0; k < NR; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 32; j++)
			{
                uint32_t *t = (uint32_t*)&matrix[k][j];
				printf("%u\t", t[i]);
			}
			printf("\n\n\n");
		}
		printf("********************************ROW**END*********************************\n");
	}

}

/**
* Print 1-d simd matrix
**/
void print1dmatrix(__m128i matrix[],int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < SLOT; j++)
		{
            uint32_t *t = (uint32_t*)&matrix[i];
			printf("%u\t",t[j]);
		}
	}
	printf("\n");
}

/**
*   bit列存储的变换必须是bit的方阵，数组的长度=每个数的bit数
**/
void invert(uint32_t *src, uint32_t *inverse)
{
	for (int i = 0; i < NB; i++)               //mask loop
	{
		for (int j = 0; j<NB; j++)         //src loop
		{
			inverse[i] += i - j>0 ? (src[j] & mask[i]) << (i - j) : (src[j] & mask[i]) >> (j - i);
		}
	}
	//print(inverse, NB);
}

/*
* Pack 4 inverse array to mach one row of simd matrix. The length of src must be 128
**/
void pack2simdrow(uint32_t *src, uint32_t inverse[SLOT][NB])
{
	uint32_t srclet[NB];
	uint32_t inverselet[NB];

	for (int i = 0; i<SLOT; i++)
	{
		memset(srclet, 0, sizeof(uint32_t)*NB); // init the srclet
		memset(inverselet, 0, sizeof(uint32_t)*NB); // init the inverselet
		memcpy(srclet, &src[i*NB], sizeof(uint32_t)*NB);
		invert(srclet, inverselet);
		//print(inverselet, NB);
		memcpy(inverse[i], inverselet, sizeof(uint32_t)*NB);
	}

}

/**
* Generating the data subject to some kind of distribution
*/
void datagenerator(uint32_t *data, int count)
{
	/*for (int i = 1; i <= count; i++)
	{
		if (i % 32 == 0 || i % 32 == 31)
			data[i - 1] = 3;
		else
			data[i - 1] = 1;
	}*/
	/*srand((unsigned)time(NULL));
	for (int i = 1; i < count; i++)
	{
		data[i - 1] = (int)(double)rand()*(double)(MAX/RAND_MAX);
		//data[i - 1] = rand();
	}*/
	for (int i = 0; i < 128; i++)
	{
		fscanf(readdata, "%d", &data[i]);
	}
}

/**
*  Read data from file
**/
/**
void datareader()
{
	readdata = fopen("/Users/zhifei/Desktop/randomnumbers_0.1billion.txt", "r+");
	int buffer[128];
}
**/

/**
* Load inverted data to a row of simd matrix
**/
void loadrow(__m128i matrix[][NB], int row, uint32_t rowdata[SLOT][NB])
{
	//print2d(rowdata,4);
	for (int i = 0; i<NB; i++)
	{
		matrix[row][i] = _mm_set_epi32(rowdata[0][i], rowdata[1][i], rowdata[2][i], rowdata[3][i]);
	}
}


/**
* Load all data into simd matrix. For now, the total number MUST BE a multiple of 128.
*/
void load2simdmatrix(__m128i matrix[][NB])
{
	
	uint32_t inverse[SLOT][NB];
	readdata=fopen("/Users/zhifei/Desktop/randomnumbers_0.1billion.txt", "r+");
    for (int i = 0; i<NR; i++) //read 128 numbers for a time
	{
		datagenerator(data+i*SLOT*NB, 128);
		pack2simdrow(data+i*SLOT*NB, inverse);
		loadrow(matrix, i, inverse);
	}
	fclose(readdata);
}


void find_init(uint32_t value)
{
	for (int i = 0; i < NR; i++)
	{
		res[i] = _mm_set_epi32(MAX, MAX, MAX, MAX);	//初始化res,让他的每一位都是1
	}

	for (int i = 0; i < NB; i++)                  //形成比较向量
	{
		int c = (value >> (31 - i)) % 2;      //获取当前比较位
		//printf("%u ", c);
		if (c == 0)
			com_value[i] = _mm_set_epi32(0, 0, 0, 0);
		else
			com_value[i] = _mm_set_epi32(MAX, MAX, MAX, MAX);      //初始化比较数组   
	}
}

void find(__m128i my[NR][32], int n, int value, int m)    //n在这里取32，value是要查找的值，m是NR
{

	__m128i tmp;									//存放临时结果
	__m128i t = _mm_set_epi32(MAX, MAX, MAX, MAX);  //辅助数组，初始化为全部都是1
    __m128i zeros = _mm_setzero_si128 ();           //0 vector
	//__m128i check = _mm_set_epi32(0, 0, 0, 0);
	//__m128i check_res;
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{ 
			tmp = _mm_andnot_si128(_mm_xor_si128(my[j][i], com_value[i]), t);  //只有当my[j][i]和Value中对应的位相同时，tmp对应的位才是1. （没有非异或函数，所以先异或，再取反）
			res[j] = _mm_and_si128(res[j], tmp); //与上一轮结果进行与操作，这样只有全部为1的最终结果才是1
			
            if (_mm_movemask_epi8(_mm_cmpeq_epi32(res[j], zeros)) == 0xffff)
            {
                //printf("BREAK\n");
                break;
            }
			

		}
	}
	//print1dmatrix(res,NR);
}


void find_without_sse(uint32_t* src, int n, int value)  //顺序查找
{

	for (int i = 0; i < n; i++)
	{
		if (src[i] == value)
		{
			res_without_sse[i] = 1;
		}
	}
}

int main()
{
	
	load2simdmatrix(matrix);

	//memset(res_without_sse, 0, sizeof(res_without_sse));      //初始化顺序查找结果数组
	//print2dmatrix(matrix);

	clock_t begin, end;
	//4294967295
	uint32_t value = 1979767798;
	//print(data,N);

	begin = clock();
	find_init(value);
	find(matrix, 32,value , NR);
	end = clock();
	printf("Use SSE：%lf \n", (double)(end - begin)/CLOCKS_PER_SEC);

	begin = clock();
	find_without_sse(data,N,value);
	end = clock();
	printf("Without SSE：%lf \n", (double)(end - begin) / CLOCKS_PER_SEC);
    
	return 0;
}