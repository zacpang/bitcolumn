#include "bitcolumn.h"
#include "print_helper.h"

// Mask of 32b and 128b
static uint32_t mask_32[32];
static __m128i mask_128[128];
static uint32_t power[32];  //The power of 2 from 1 to 32. The pow[32] is 2^32

__m128i res[N_ROWS];								//保存运用SSE时，比较的结果
__m128i matrix[N_ROWS][N_BITS];
__m128i com_value[N_BITS];                         //由将要比较的数得到，生成比较的数组
uint32_t data[N];
FILE *readdata;

uint32_t res_set[10000];

/*
 * Init the environment. The 32b-mask, 128b-mask and the power array of 2
 from 1 to 31.
*/
void init()
{
    uint32_t powof2 = _2POW31_;
    for (int j = 0; j < 32; j++)
    {
        mask_128[j] = _mm_set_epi32(powof2, 0, 0, 0);
        mask_128[32 + j] = _mm_set_epi32(0, powof2, 0, 0);
        mask_128[64 + j] = _mm_set_epi32(0, 0, powof2, 0);
        mask_128[96 + j] = _mm_set_epi32(0, 0, 0, powof2);

        mask_32[j] = powof2;
        power[31-j] = powof2;
        powof2 = powof2>>1;
    }

}

/*
 * If is the v is full of 0s, return 1; 0 otherwise
 */
int is_zero(__m128i* v)
{
    static __m128i zeros = _mm_setzero_si128();
    return _mm_movemask_epi8(_mm_cmpeq_epi32(*v, zeros)) == 0xffff;
}

/*
*    Invert the array to bit-column-way. The length of src and inverse must be 32!
*
*/
void uint32_invert(uint32_t *src, uint32_t *inverse)
{
	for (int i = 0; i < N_BITS; i++)               //mask loop
	{
		for (int j = 0; j<N_BITS; j++)         //src loop
		{
			inverse[i] += i - j>0 ? (src[j] & mask_32[i]) << (i - j) : (src[j] & mask_32[i]) >> (j - i);
		}
	}
	//print(inverse, N_BITS);
}


/*
* Pack 4 inverse array to mach one row of simd matrix. The length of src must be 128
*/
void pack2simdrow(uint32_t *src, uint32_t inverse[][N_BITS])
{
	uint32_t srclet[N_BITS];
	uint32_t inverselet[N_BITS];

	for (int i = 0; i<N_SLOTS; i++) //For 32bit number, N_SLOTS would be 4;
	{
		memset(srclet, 0, sizeof(uint32_t)*N_BITS); // init the srclet
		memset(inverselet, 0, sizeof(uint32_t)*N_BITS); // init the inverselet
		memcpy(srclet, &src[i*N_BITS], sizeof(uint32_t)*N_BITS);
		invert(srclet, inverselet);
		//print(inverselet, N_BITS);
		memcpy(inverse[i], inverselet, sizeof(uint32_t)*N_BITS);
	}

}

/**
* Generating the data subject to some kind of distribution
*/
void datagenerator(uint32_t *data, int count)
{
	for (int i = 0; i < count; i++)
	{
		data[i] = i;
	}
	/*srand((unsigned)time(NULL));
	for (int i = 1; i < count; i++)
	{
	data[i - 1] = (int)(double)rand()*(double)(MAX/RAND_MAX);
	//data[i - 1] = rand();
	}*/
	/*for (int i = 0; i < 128; i++)
	{
		fscanf(readdata, "%d", &data[i]);
	}*/
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

void loadrow(__m128i matrix[][N_BITS], int row, uint32_t rowdata[N_SLOTS][N_BITS])
{
	//print2d(rowdata,4);
	for (int i = 0; i<N_BITS; i++)
	{
		matrix[row][i] = _mm_set_epi32(rowdata[0][i], rowdata[1][i], rowdata[2][i], rowdata[3][i]);
	}
}


/**
* Load all data into simd matrix. For now, the total number of data MUST BE a multiple of 128.
*/
void load2simdmatrix(__m128i matrix[][N_BITS])
{

	uint32_t inverse[N_SLOTS][N_BITS];
	readdata = fopen("/Users/zhifei/Desktop/uniform_0.1billion.txt", "r+");
	for (int i = 0; i<N_ROWS; i++) //read 128 numbers for a time
	{
		datagenerator(data + i*N_SLOTS*N_BITS, 128);
		pack2simdrow(data + i*N_SLOTS*N_BITS, inverse);
		loadrow(matrix, i, inverse);
	}
	fclose(readdata);
}

/**
 * 1. Convert the value to simd array. Range will be converted into a simd array whose length is # of bits
 * shared in common.
 * 2. A 128-bit mask will be crated.
 * 3. A power array which contains the power of 2 from 1 to 31 for the combination from bits to integer.
 *
 **/
int find_init(uint32_t left, uint32_t right)
{
	int count = 0;
	for (int i = 0; i < N_ROWS; i++)
	{
		res[i] = _mm_set_epi32(MAX, MAX, MAX, MAX);	//初始化res,让他的每一位都是1
	}

	for (int i = 0; i < N_BITS; i++)                  //形成比较向量
	{
		int c = (left >> (31 - i)) % 2;      //获取当前比较位
		int d = (right >> (31 - i)) % 2;
		//printf("%u ", c);
		if (c == d)
		{
			count++; //The # of bits shared in the common prefix of left and right
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

	return count;
}

int res_index = 0;
void check(int r, int com_res[], int num,int left,int right) //检查数据是否在left和right之间
{
	int pri_value = 0;   //记录要还原的数据
	for (int i = 0; i < num; i++)
	{
		pri_value = 0;
		for (int j = 0; j < N_BITS; j++)
		{
			__m128i tmp = _mm_and_si128(matrix[r][j], mask_128[com_res[i]]);   //如果matrix[r][j]中第com_res[i]位是1，tmp中对应位也是1
            uint64_t *t = (uint64_t*)&tmp;
			if (t[0] != 0 || t[1] != 0)  //判断tmp是否是全0
			{
				pri_value += power[31 - j];
			}
		}
		//printf("pri_value: %u \n", pri_value);
		if (pri_value >= left && pri_value <= right)
		{
			res_set[res_index] = pri_value;  //res_set存放最终结果
			res_index++;
		}
	}
}
/**
 * n is number of bit
 */
void find(__m128i my[N_ROWS][32], int n, int left,int right, int m)    //n在这里取32，value是要查找的值，m是N_ROWS
{

	__m128i tmp;									      //存放临时结果
	__m128i t = _mm_set_epi32(MAX, MAX, MAX, MAX);        //辅助数组，初始化为全部都是1
	__m128i zeros = _mm_setzero_si128();                  //0 vector
	int count = find_init(left, right);  // count is the # of bits shared in common;
	for (int j = 0; j < m; j++)  //m is the number of rows of data matrix
	{
		for (int i = 0; i < count; i++)
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
    /**
     * Check the mask_128. If mask_128[i] != 0, check this.
     **/
	int com_res[128];    //存放待进一步检查的位
	int next = 0;
	for (int i = 0; i < N_ROWS; i++)
	{
		for (int j = 0; j < N_BITS*N_SLOTS; j++)
		{
			__m128i tmp = _mm_and_si128(res[i], mask_128[j]);
            uint64_t *t = (uint64_t*)&tmp;
			if (t[0] != 0 || t[1] != 0)
			{
				com_res[next] = j;        //记录要进一步检查的位
				next++;
			}
		}
        //i 代表检查到第几个向量，next表示com_res的长度。com_res记录这个向量不为0的位
		check(i, com_res, next,left,right);
		/*printf("\n");
		for (int k = 0; k < next; k++)
		{
			printf("%d - ",com_res[k]);
		}
		printf("\n");*/
	}
	//print1dmatrix(res,N_ROWS);
	/*for (int i = 0; i < N_ROWS; i++)
	{
		for (int j = 0; j < N_SLOTS; j++)
		{
			if (res[i].m128i_i32[j] == 0)
			{
				continue;
			}
			else
			{
				//printf("\n == %u ==\n", res[i].m128i_i32[j]);
				for (int k = 0; k < N_BITS; k++)
				{
					int c = (res[i].m128i_i32[j] >> k) % 2;
					if (c == 0)
						continue;
					else
					{
						//printf("%d %d %d %d\n ", i, j, k, data[i * 128 + (N_SLOTS - 1 - j)*N_BITS + (N_BITS - 1 - k)]);
						//if (data[i * 128 + (N_SLOTS - 1 - j)*N_BITS + (N_BITS - 1 - k)] < left || data[i * 128 + (N_SLOTS - 1 - j)*N_BITS + (N_BITS - 1 - k)] > right)
						int pri_value = 0;
						for (int l = 0; l < N_BITS; l++)
						{
							int flag= (uint32_t)(matrix[i][l].m128i_i32[j] >> k) % 2 ;
							//printf("%d - %u %d %d %d \n", flag, matrix[i][l].m128i_i32[j],i,l,j);
							pri_value +=flag * power[31 - l];
						}
						//printf("\npri_value : %d \n",pri_value);
						if (pri_value < left || pri_value > right )
						{
							//printf("\n pri_value: %d    left: %d     right:%d sub\n", pri_value, left,right);
							res[i].m128i_i32[j] = res[i].m128i_i32[j] - power[k];
						}
						else
						{
							res_num[count_res] = pri_value;
							count_res++;
						}
					}
				}
			}
		}
	}*/
}


void find_without_sse(uint32_t* src, int n, int left,int right)  //顺序查找
{

	for (int i = 0; i < n; i++)
	{
		if (src[i] >= left && src[i] <= right)
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
	uint32_t left = 47;
	uint32_t right = 90;
	//print(data,N);

	//print2dmatrix(matrix);

	begin = clock();

	find(matrix, 32, left,right, N_ROWS);
	end = clock();
	printf("Use SSE：%lf \n", (double)(end - begin) / CLOCKS_PER_SEC);
	//print1dmatrix(res,N_ROWS);
	printf("final number: %d\n",res_index);
	for (int j = 0; j < res_index; j++)
	{
		printf("%u ", res_set[j]);
	}
	begin = clock();
	find_without_sse(data, N, left,right);
	end = clock();
	printf("Without SSE：%lf \n", (double)(end - begin) / CLOCKS_PER_SEC);
	//print1dmatrix(mask_128, 128);
	return 0;
}
