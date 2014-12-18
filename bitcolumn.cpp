#include "bitcolumn.h"
#include "print_helper.h"

static uint32_t mask_32[32];
static __m128i mask_128[128];
static uint32_t power[32];  //The power of 2 from 1 to 32. The pow[32] is 2^32
static __m128i res_mask[N_ROWS];	//保存运用SSE时，比较的结果
static __m128i ones, zeros;  //A simd vector with all 1s or 0s

static __m128i matrix[N_ROWS][N_BITS];
uint32_t res_set[10000]; //存放最终结果的10进制的数字，将来可以使用动态链表存放，为方便起见，开一个10^4数组
int res_index = 0;

/*
 * Init the environment.
 0) init the 1s and 0s vector
 1) Generate a 32 bit mask, 128b mask and the power array of 2 from 1 to 31;
 2) Init the result mask;
 */
void init()
{
    ones = _mm_set_epi32(UINT32MAX, UINT32MAX, UINT32MAX, UINT32MAX);
	zeros = _mm_setzero_si128();

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

    for(int i=0; i<N_ROWS; i++)
    {
        res_mask[i] = ones;
    }
}

/*
 * If is the v is full of 0s, return 1; 0 otherwise ?
 */
int is_zero(__m128i* v)
{
    return _mm_movemask_epi8(_mm_cmpeq_epi32(*v, zeros)) == 0xffff;
}

/*
*  Invert the array to bit-column-way. The length of src and inverse MUST BE 32!
*
*/
void uint32_invert(uint32_t *src, uint32_t *inverse)
{
	for (int i = 0; i < N_BITS; i++)               //mask loop
	{
		for (int j = 0; j<DATATYPE_LEN; j++)         //src loop
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

	uint32_t srclet[DATATYPE_LEN];
	uint32_t inverselet[N_BITS];

	for (int i = 0; i<N_SLOTS; i++) //For 32bit number, N_SLOTS would be 4;
	{
		memset(srclet, 0, sizeof(uint32_t)*DATATYPE_LEN); // init the srclet
		memset(inverselet, 0, sizeof(uint32_t)*N_BITS); // init the inverselet
		memcpy(srclet, &src[i*DATATYPE_LEN], sizeof(uint32_t)*DATATYPE_LEN);
		uint32_invert(srclet, inverselet);
		memcpy(inverse[i], inverselet, sizeof(uint32_t)*N_BITS);
	}
}

/*
 * Fill one row of simd matrix with inverse array.
 */
void loadrow(__m128i matrix[][N_BITS], int row, uint32_t rowdata[N_SLOTS][N_BITS])
{
	for (int i = 0; i<N_BITS; i++)
	{
		matrix[row][i] = _mm_set_epi32(rowdata[3][i], rowdata[2][i], rowdata[1][i], rowdata[0][i]);
	}
}


/**
* Load all data into simd matrix. For now, the total number of data MUST BE a multiple of 128.
*/
void load2simdmatrix(__m128i matrix[][N_BITS])
{
	uint32_t inverse[N_SLOTS][N_BITS];
    uint32_t buffer[V_LEN];
    int buffer_idx = 0;
    uint32_t row_counter = 0;
    FILE *datafile = fopen("/Users/zhifei/Desktop/test.txt", "r+");

    for(int i=0; i<N; i++)
    {
        fscanf(datafile,"%u",&buffer[buffer_idx++]);
        if(buffer_idx == V_LEN)
        {
            pack2simdrow(buffer, inverse);
            loadrow(matrix, row_counter++, inverse);
            buffer_idx = 0;
        }
    }
}

/**
 Convert single value into bit-column-store way. The length of v should be the # of bits of value
 */
void convert_bitcolumn_single(uint32_t value, __m128i *v)
{
    uint32_t set_value;

    for(int i=0; i<N_BITS; i++)
    {
        set_value = ((value >> (31 - i)) % 2) * UINT32MAX;
        v[i] = _mm_set_epi32(set_value, set_value, set_value, set_value);
    }
}

/**
 * 将两个数的共同高位bit串，转换成SIMD数组，返回真实的共同bit个数.如果共同bit数大于length,返回length。
 * （注：这样做是因为，以后可能uint32_t的数只存储16位到SIMD，所以如果共同bit数不能大于N_BITS）
 */
int convert_bitcolumn_double(uint32_t left, uint32_t right, __m128i *v)
{
    int real_length = 0;
    int flag1, flag2;
    uint32_t set_value;
    for(int i=0; i<N_BITS; i++)
    {
        flag1 = (left >> (31 - i)) % 2;
        flag2 = (right >> (31 - i)) % 2;
        if(flag1 != flag2)
        {
            break;

        } else {
            set_value = flag2 * UINT32MAX;
            v[i] = _mm_set_epi32(set_value, set_value, set_value, set_value);
            real_length++;
        }
    }

    return real_length;
}

/*
 tofind is the _m128i array which is converted by the value to find. The length
 indicates the actual length of SIMD array. 单点查询和范围查询均调用此函数
 */
void do_search(__m128i simd_mtx[N_ROWS][N_BITS], __m128i *to_find, int prefix_length)
{
    __m128i tmp;
    for(int row=0; row<N_ROWS; row++)
    {
        for(int i=0; i<prefix_length; i++)
        {
            tmp = _mm_andnot_si128(_mm_xor_si128(simd_mtx[row][i], to_find[i]), ones); //只有当simd_mtx[row][i]和tofind[i]中对应的位相同时，tmp对应的位才是1.
            res_mask[row] = _mm_and_si128(res_mask[row], tmp);

            if (is_zero(&res_mask[row]))
			{
				break;
			}
        }
    }
}

/**
* 单点查询
*/
void single_search(uint32_t value)
{
    __m128i *to_find = (__m128i *)malloc(sizeof(__m128i)*N_BITS);
    convert_bitcolumn_single(value, to_find);
    do_search(matrix, to_find, N_BITS); //执行完,res_mask中即得到对应元素的掩码。
}

/**
* 范围查询
*/
void range_search(uint32_t left, uint32_t right)
{
    __m128i *to_find = (__m128i *)malloc(sizeof(__m128i)*N_BITS);

    int real_length = convert_bitcolumn_double(left, right, to_find);

    do_search(matrix, to_find, real_length); //执行完,res_mask中即得到对应元素的掩码。
    
    
    /**
     * 如果共同前缀长度小于8，拼接
     */
    
}

/**
*将位存储的数字组装成十进制的数字，然后和范围短点比较。
*/
void check(int r, int com_res[], int num,int left,int right) //检查数据是否在left和right之间
{
	int pri_value = 0;   //记录要还原的数据
	uint64_t *t;
	for (int i = 0; i < num; i++) //num is the length of com_res
	{
		pri_value = 0;
		for (int j = 0; j < N_BITS; j++)
		{
			__m128i tmp = _mm_and_si128(matrix[r][j], mask_128[com_res[i]]);   //如果matrix[r][j]中第com_res[i]位是1，tmp中对应位也是1
            t = (uint64_t*)&tmp;
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
*　此函数用于数据全部存储于simd_matrix的情况
*/
void range_validate_v1(uint32_t left, uint32_t right)
{
    int res_validate_bit_loc[V_LEN];    //存放待进一步检查的位
	int rvbl_idx = 0;
    int flag;
    __m128i tmp;
	for (int row = 0; row < N_ROWS; row++)
	{
		for (int j = 0; j <V_LEN; j++)
		{
            tmp = _mm_and_si128(res_mask[row], mask_128[j]);
            flag = is_zero(&tmp);
            res_validate_bit_loc[rvbl_idx] += j*flag;
            rvbl_idx += flag;
		}
        check(row, res_validate_bit_loc, rvbl_idx, left, right);
	}
}

void clear_res_set()
{
    memset(res_set, 0, sizeof(uint32_t)*10000);
    res_index = 0;
}




int main()
{
    init(); //This cant be deleted
    
	load2simdmatrix(matrix);
    
    clear_res_set();
    
    range_search(3, 5);
    
    
    //print2dmatrix(matrix);

	//memset(res_without_sse, 0, sizeof(res_without_sse));      //初始化顺序查找结果数组
	//print2dmatrix(matrix);
/*
	clock_t begin, end;

	uint32_t left = 47;
	uint32_t right = 90;

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
	*/
	return 0;
}
