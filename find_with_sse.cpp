#include <stdio.h>
#include <inttypes.h>
#include <memory.h>
#include <xmmintrin.h>
#include <time.h>
#include <stdlib.h>
#include <emmintrin.h>

#define NB 32                                  //number of bits of uint32_t
#define SLOT (128/NB)
#define MAX 4294967295                         //the MAX value of uint32_t
#define N 1280000                   //the # of number
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

__m128i com_mask[128];


__m128i res[NR];								//保存运用SSE时，比较的结果
uint32_t res_without_sse[N];                    //保存顺序查询时，比较的结果
__m128i matrix[NR][NB];
__m128i com_value[NB];                         //由value得到的，用于比较的数组
uint32_t data[N];
FILE *readdata;
uint32_t power[32];
uint32_t res_set[100000];
/**
*  Pirnt 1-d array
**/
void print(uint32_t *a, int count)
{
	int i;
	for (i = 0; i<count; i++)
	{
		if (i % 32 == 0)
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
void print1dmatrix(__m128i matrix[], int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < SLOT; j++)
		{
			uint32_t *t = (uint32_t*)&matrix[i];
			printf("%u\t", t[j]);
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
	/*for (int i = 0; i < count; i++)
	{
		data[i] = i;
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
	readdata = fopen("D:\\文档\\data\\randomnumbers_0.1billion.txt", "r+");
	for (int i = 0; i<NR; i++)   //read 128 numbers for a time
	{
		datagenerator(data + i*SLOT*NB, 128);
		pack2simdrow(data + i*SLOT*NB, inverse);
		loadrow(matrix, i, inverse);
	}
	fclose(readdata);
}


int find_init(uint32_t left, uint32_t right)
{
	int count = 0;
	for (int i = 0; i < NR; i++)
	{
		res[i] = _mm_set_epi32(MAX, MAX, MAX, MAX);	//初始化res,让他的每一位都是1
	}

	for (int i = 0; i < NB; i++)                  //形成比较向量
	{
		int c = (left >> (31 - i)) % 2;      //获取当前比较位
		int d = (right >> (31 - i)) % 2;
		if (c == d)
		{
			count++;
			if (c == 0)
				com_value[i] = _mm_set_epi32(0, 0, 0, 0);
			else
			{
				com_value[i] = _mm_set_epi32(MAX, MAX, MAX, MAX);      //初始化比较数组   
			}
		}
		else
			break;
	}
	//printf("\n %d \n",count);
	uint32_t sum = 1;

	for (int i = 0; i < 4; i++)   //初始化com_mask
	{
		sum = 2147483648;
		for (int j = 0; j < 32; j++)
		{
			switch (i)
			{
			case 0:com_mask[i * 32 + j] = _mm_set_epi32(sum, 0, 0, 0); break;
			case 1:com_mask[i * 32 + j] = _mm_set_epi32(0, sum, 0, 0); break;
			case 2:com_mask[i * 32 + j] = _mm_set_epi32(0, 0, sum, 0); break;
			case 3:com_mask[i * 32 + j] = _mm_set_epi32(0, 0, 0, sum); break;
			default:
				break;
			}
			sum /= 2;
		}
	}

	sum = 1;
	for (int i = 0; i < 32; i++)  //初始化power数组
	{
		power[i] = sum;
		sum *= 2;
	}
	return count;
}

int count_withsse = 0;
void check(int r, int com_res[], int num, int left, int right) //检查数据是否在left和right之间
{
	int pri_value = 0;   //记录要还原的数据
	for (int i = 0; i < num; i++)
	{
		pri_value = 0;
		for (int j = 0; j < NB; j++)
		{
			__m128i tmp = _mm_and_si128(matrix[r][j], com_mask[com_res[i]]);   //如果matrix[r][j]中第com_res[i]位是1，tmp中对应位也是1
			if (tmp.m128i_i32[0] != 0 || tmp.m128i_i32[1] != 0 || tmp.m128i_i32[2] != 0 || tmp.m128i_i32[3] != 0)  //判断tmp是否是全0
			{
				pri_value += power[31 - j];
			}
		}
		//printf("pri:%d\n",pri_value);
		if (pri_value >= left && pri_value <= right)
		{
			res_set[count_withsse] = pri_value;  //res_set存放最终结果       
			count_withsse++;
		}
	}
}

void find(__m128i my[NR][32], int n, int left, int right, int m)    //n在这里取32，value是要查找的值，m是NR
{

	__m128i tmp;									      //存放临时结果
	__m128i t = _mm_set_epi32(MAX, MAX, MAX, MAX);        //辅助数组，初始化为全部都是1
	__m128i zeros = _mm_setzero_si128();                  //0 vector
	int count = find_init(left, right);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < count; i++)
		{
			tmp = _mm_andnot_si128(_mm_xor_si128(my[j][i], com_value[i]), t);  //只有当my[j][i]和Value中对应的位相同时，tmp对应的位才是1. （没有非异或函数，所以先异或，再取反）
			res[j] = _mm_and_si128(res[j], tmp); //与上一轮结果进行与操作，这样只有全部为1的最终结果才是1

			if (_mm_movemask_epi8(_mm_cmpeq_epi32(res[j], zeros)) == 0xffff)
			{
				break;
			}
		}
	}
	int com_res[128];    //存放待进一步检查的位
	int next = 0;        //向com_res中存数时候的下标
	int pri_value = 0;
	for (int i = 0; i < NR; i++)
	{
		next = 0;
		for (int j = 0; j < NB*SLOT; j++)
		{
			__m128i tmp = _mm_and_si128(res[i], com_mask[j]);
			if (tmp.m128i_i64[0] != 0 || tmp.m128i_i64[1] != 0)
			{
				com_res[next] = j;        //记录要进一步检查的位
				next++;
			}
		}

		check(i, com_res, next, left, right);
		/*printf("\n");
		for (int k = 0; k < next; k++)
		{
		printf("%d - ",com_res[k]);
		}
		printf("\n");*/
	}
}

/**
* find the data without see,just check the data from left to right
**/
int count_without = 0;   //record the number of nodes that meet the request
void find_without_sse(uint32_t* src, int n, int left, int right)    //顺序查找
{

	for (int i = 0; i < n; i++)
	{
		if (src[i] >= left && src[i] <= right)
		{
			res_without_sse[i] = 1;
			count_without++;
		}
	}
}



int main()
{

	load2simdmatrix(matrix);    //将数据转移到matrix中

	//print2dmatrix(matrix);

	clock_t begin, end;
	uint32_t left = 0;
	uint32_t right = 13000;
	//print(data,N);

	begin = clock();
	find(matrix, 32, left, right, NR);
	end = clock();
	printf("Use SSE：%lf \n", (double)(end - begin) / CLOCKS_PER_SEC);
	//print1dmatrix(res,NR);
	printf("final number: %d\n", count_withsse);   //打印出符合要求的节点数目

	begin = clock();
	find_without_sse(data, N, left, right);
	end = clock();
	printf("Without SSE：%lf \n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("final number: %d\n", count_without);  //打印出顺序查找时符合要求的节点数目

	return 0;
}