#include <stdio.h>
#include <inttypes.h>
#include <memory.h>
#include <emmintrin.h>
#include <time.h>


#define NB 32                                  //number of bits of uint32_t
#define SLOT (128/NB)
#define MAX 4294967295                         //the MAX value of uint32_t
#define N 128                                  //the # of number
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


__m128i my[NR][32];                              //保存重构之后插入的数字
uint32_t src[N];								//初始插入数组
__m128i res[NR];								//保存运用SSE时，比较的结果
uint32_t res_without_sse[N];                    //保存顺序查询时，比较的结果

/**
*bit列存储的变换必须是bit的方阵，数组的长度=每个数的bit数
**/
void invert(uint32_t *src, uint32_t *inverse, int count)
{
	if (count != NB)
	{
		printf("The length of the array must be equal to number of bits");
        return;
	}
	
	memset(inverse, 0, sizeof(uint32_t)*count); // init the inverse
    
    int i, j;
    
    //test
    
    for (i = 0; i < NB; i++)
    {
        printf("%u\t", src[i]);
    }
    printf("\n\n\n");
    
    //test end
    
    
    
    
	
	
	for (i = 0; i < NB; i++)               //mask loop
	{
		for (j = 0; j<count; j++)         //src loop
		{
			inverse[i] += i-j>0 ? (src[j] & mask[i]) << (i-j) : (src[j] & mask[i]) >> (j-i);
		}
	}
    
    //test
    /**
    for (i = 0; i < NB; i++)
    {
        printf("%u\t", inverse[i]);
    }
    printf("\n\n\n");**/
    
    //test end
}

/**
*read 128 numbers and seal into uint32_t[4][32] to fill 1 row of _m128i[X][32]
*The length of the src must be 128
pack的每一列填充一个_m128i向量，一共32列，可以填充32个向量。在向量的二维数组里，可以填充一行。
**/
/*void pack2row(uint32_t *src, uint32_t pack[][4][NB])
{
	int i = 0;
	int j = 0;
	int loop = 4; //将src数组
	uint32_t srclet[NB];
	uint32_t inverse[NB];
	for (j = 0; j < NR; j++)
	{
		for (i = 0; i < loop; i++)
		{
			memset(srclet, 0, sizeof(srclet));
			memset(inverse, 0, sizeof(inverse));
			memcpy(srclet, &src[(j*4+i)*NB], sizeof(srclet));
			invert(srclet, inverse, NB);
			memcpy(&pack[j][i], inverse, sizeof(inverse));
		}
	}
}*/


void pack2row(uint32_t *src, __m128i my[NR][32])
{
    int i, j, k;
	
	uint32_t srclet[NB];
	uint32_t inverse[NB];
	for (j = 0; j < NR; j++)
	{
		for (i = 0; i < SLOT; i++)
		{
			memset(srclet, 0, sizeof(srclet));
			memset(inverse, 0, sizeof(inverse));
			memcpy(srclet, &src[(j * 4 + i)*NB], sizeof(srclet));
			invert(srclet, inverse, NB);
			for (k = 0; k < NB; k++)
			{
				//memcpy(&my[j][k].m128i_i32[i], inverse + k, sizeof(uint32_t));    //把转换后的数字写入my中
                memcpy(&my[j][k]+NB*i, inverse+k, sizeof(uint32_t));    //把转换后的数字写入my中
                //my[j][k].m128i_i32[i] = inverse[k]
            }
		}
	}
}



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
	/*printf("\n");
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			printf("%u  ", res[j].m128i_u32[i]);
		}
	}*/
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
	
	/**
	uint32_t src[NB] = {1,1,1,1,1,1,1,1,
	1,1,1,1,1,1,1,1,
	1,1,1,1,1,1,1,1,
	1,1,1,1,1,1,1,7};
	uint32_t inverse[NB];
	memset(inverse,0,sizeof(inverse));

	invert(src,inverse,NB);


	int i = 0;
	for(i=0; i<NB; i++)
	{
	printf("%ld\n",inverse[i]);
	}
	**/

	 //读入的128个数
	//uint32_t pack[(N + 127) / 128][4][NB];
	int i, j, k;

	for (i = 1; i <= N; i++) //模拟128个数字
	{
        src[i-1] = 1;
    }
	
	
	memset(res_without_sse, 0, sizeof(res_without_sse));      //初始化顺序查找结果数组
	pack2row(src, my);
    
    
    
	//int loop = (N+31)/32;
	//输出128个数
	
	/**
	for (k = 0; k < NR; k++)
	{
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 32; j++)
			{
                __m128i tmp = my[k][j];
                
				printf("%lld\t", tmp[i]);
			}
			printf("\n\n\n");
		}
		printf("__________________\n");
	}**/
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
