#include <stdio.h>
#include <inttypes.h>
#include <memory.h>
#include <emmintrin.h>
#include <time.h>


#define NB 32                                  //number of bits of uint32_t
#define SLOT (128/NB)
#define MAX 4294967295                         //the MAX value of uint32_t
#define N 256                                  //the # of number
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


/**
 *  Pirnt 1-d array
 **/
void print(uint32_t *a, int count)
{
    int i;
    for (i=0; i<count; i++)
    {
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
    for (i=0; i<d; i++)
    {
        for (j=0; j<32; j++)
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
void printmatrix(__m128i matrix[][NB])
{
   
    for (int k = 0; k < NR; k++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 32; j++)
            {
                __m128i tmp = matrix[k][j];
                printf("%u\t", tmp[i]);
            }
            printf("\n\n\n");
        }
        printf("********************************ROW**END*********************************\n");
    }
    
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
			inverse[i] += i-j>0 ? (src[j] & mask[i]) << (i-j) : (src[j] & mask[i]) >> (j-i);
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
    
    for(int i = 0; i<SLOT; i++)
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
    for (int i=1; i<=count; i++)
    {
        if(i%32 == 0)
            data[i-1] = 3;
        else
            data[i-1] = 1;
    }
}

/**
 *  Read data from file
 **/
void datareader(){}


/**
 * Load inverted data to a row of simd matrix
 **/
void loadrow(__m128i matrix[][NB], int row, uint32_t rowdata[SLOT][NB])
{
    //print2d(rowdata,4);
    for(int i=0; i<NB; i++)
    {
        matrix[row][i] = _mm_set_epi32(rowdata[0][i], rowdata[1][i], rowdata[2][i], rowdata[3][i]);
    }
}


/**
 * Load all data into simd matrix. For now, the total number MUST BE a multiple of 128.
 */
void load2simdmatrix(__m128i matrix[][NB])
{
    uint32_t data[N];
    uint32_t inverse[SLOT][NB];
    
    for (int i=0; i<NR; i++) //read 128 numbers for a time
    {
        datagenerator(data, 128);
        pack2simdrow(data, inverse);
        loadrow(matrix, i, inverse);
    }
}


/**
void pack2row(uint32_t *src, __m128i my[NR][32])
{
    int i, j, k;

	uint32_t srclet[NB];
	uint32_t inverse[NB];
	for (j = 0; j < NR; j++)
	{
		for (i = 0; i < SLOT; i++)
		{
			memset(srclet, 0, sizeof(uint32_t)*NB);
			memset(inverse, 0, sizeof(uint32_t)*NB);
			memcpy(srclet, &src[(j*4+i)*NB], sizeof(uint32_t)*NB);
            
			invert(srclet, inverse);
			for (k = 0; k < NB; k++)
			{
				memcpy(&my[j][k].m128i_i32[i], inverse + k, sizeof(uint32_t));    //把转换后的数字写入my中
                //memcpy(&my[j][k]+NB*i, inverse+k, sizeof(uint32_t));    //把转换后的数字写入my中
                //my[j][k].m128i_i32[i] = inverse[k]
            }
		}
	}
}
**/


void find(__m128i my[NR][32], int n, int value, int m)    //n在这里取32，value是要查找的值，m是NR
{
	
	__m128i Value;        //由value得到的，用于比较的数组
	for (int i = 0; i < NR; i++)
	{
		res[i] = _mm_set_epi32(MAX, MAX, MAX, MAX);	//初始化res,让他的每一位都是1
	}
	__m128i tmp;									//存放临时结果				
	__m128i t = _mm_set_epi32(MAX, MAX, MAX, MAX);  //辅助数组，初始化为全部都是1
	//printf("%u ===",res.m128i_u32[0]);
	
	for (int j = 0; j < m;j++)
	{
		for (int i = 0; i < n; i++)
		{
			int c = (value >> (31 - i)) % 2;      //获取当前比较位
			//printf("%u ", c);
			if (c == 0)
				Value = _mm_set_epi32(0, 0, 0, 0);
			else
				Value = _mm_set_epi32(MAX, MAX, MAX, MAX);      //初始化比较数组   
			tmp = _mm_andnot_si128(_mm_xor_si128(my[j][i], Value), t);  //只有当my[j][i]和Value中对应的位相同时，tmp对应的位才是1. （没有非异或函数，所以先异或，再取反）
			res[j] = _mm_and_si128(res[j], tmp); //与上一轮结果进行与操作，这样只有全部为1的最终结果才是1
		}
	}
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
    __m128i matrix[NR][NB];
    
    load2simdmatrix(matrix);
    
	//memset(res_without_sse, 0, sizeof(res_without_sse));      //初始化顺序查找结果数组
	
    printmatrix(matrix);

	
    /**
	clock_t begin, end;
	begin = clock();
    find(my,32,7,NR);
	end = clock();
	printf("Use SSE：%lf \n", (double)(end - begin));

	begin = clock();
	find_without_sse(src,N,7);
	end = clock();
	printf("Without SSE：%lf \n", (double)(end - begin));**/
	return 0;
}
