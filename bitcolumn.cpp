#include "bitcolumn.h"
#include "print_helper.h"

static uint32_t mask_32[32];
static __m128i mask_128[128];
static uint32_t mask_groupby8[4]; //8位一组，分别为1，共4个数字
static uint32_t mask_groupby16[2]; //16位一组，分别为1，共2个数字
static uint32_t row_mask[4]; //用于得到后面BITS_IN_ARY位对应的数字
static uint32_t power[32];  //The power of 2 from 1 to 32. The pow[32] is 2^32
static __m128i res_mask[N_ROWS];	//保存运用SSE时，比较的结果
static __m128i ones, zeros;  //A simd vector with all 1s or 0s
static uint8_t array_data[N][CACHE_LINE_SIZE];
static __m128i matrix[N_ATTR][N_ROWS][N_BITS];
//static __m128i matrix[N_ATTR][N_ROWS][N_BITS];


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
    
    mask_groupby16[0] = 0xffff0000;
    mask_groupby16[1] = 0x0000ffff;
    
    mask_groupby8[0] = 0xff000000;
    mask_groupby8[1] = 0x00ff0000;
    mask_groupby8[2] = 0x0000ff00;
    mask_groupby8[3] = 0x000000ff;


	row_mask[0] = 0x000000ff;
	row_mask[1] = 0x0000ffff;
	row_mask[2] = 0x00ffffff;
	row_mask[3] = 0xffffffff;
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
void uint32_invert(uint32_t *src, uint32_t *inverse, int begin_index)
{
	for (int i = 0; i < N_BITS; i++)               //mask loop
	{
		for (int j = 0; j<DATATYPE_LEN; j++)         //src loop
		{
			inverse[i] += i - j>0 ? (src[j] & mask_32[i]) << (i - j) : (src[j] & mask_32[i]) >> (j - i);
		}
	}
	
	
	print(inverse, N_BITS);
}


/*
* Pack 4 inverse array to mach one row of simd matrix. The length of src must be V_LEN 128
*/
void pack2simdrow(__m128i matrix[N_ATTR][N_ROWS][N_BITS], int row, uint32_t src[V_LEN][N_ATTR])
{
	
	uint32_t src_with_attr[V_LEN];
	uint32_t inverse[N_SLOTS][N_BITS];
	uint32_t srclet[DATATYPE_LEN];
	uint32_t inverselet[N_BITS];
	
	for (int i = 0; i < N_ATTR; i++)
	{
		for (int j = 0; j < V_LEN; j++)
		{
			src_with_attr[j] = src[i][j];
		} //Get the src with specific attribute
			
		for (int m = 0; m<N_SLOTS; m++) //For 32bit number, N_SLOTS would be 4;
		{
			memset(srclet, 0, sizeof(uint32_t)*DATATYPE_LEN); // init the srclet
			memset(inverselet, 0, sizeof(uint32_t)*N_BITS); // init the inverselet
			memcpy(srclet, &src_with_attr[m*DATATYPE_LEN], sizeof(uint32_t)*DATATYPE_LEN);
			uint32_invert(srclet, inverselet, row*V_LEN + m*DATATYPE_LEN);
			memcpy(inverse[m], inverselet, sizeof(uint32_t)*N_BITS);
		}

		for (int k = 0; k<N_BITS; k++)
		{
			matrix[i][row][k] = _mm_set_epi32(inverse[3][k], inverse[2][k], inverse[1][k], inverse[0][k]);
		}

	}
    
}

/**
* Load all data into simd matrix. For now, the total number of data MUST BE a multiple of 128.
*/
void laod_all_data(__m128i matrix[][N_ROWS][N_BITS])
{
	
	uint32_t buffer[V_LEN][N_ATTR];
	uint32_t buffer_idx = 0;
    uint32_t row_counter = 0;
	uint32_t record[N_ATTR];
	memset(record, 0, sizeof(uint32_t)* N_ATTR);
    FILE *datafile = fopen("mutiple_attributes.txt", "r+");
    
    for(int i = 0; i < N; i++) //process 128 data once a time
    {
		buffer_idx++;
        if(buffer_idx == V_LEN)
        {
			for (int i = 0; i<V_LEN; i++)
			{
				for (int j = 0; j<N_ATTR; j++)
				{
					printf("%u\t", buffer[i][j]);
				}
				printf("\n");
			}
            pack2simdrow(matrix, row_counter++, buffer);
            buffer_idx = 0;
        }
        


		for (int m = 0; m < N_ATTR; m++)
		{
			for (int j = 0; j < N_BYTE_INARRAY; j++) //将所有数据处理存储到array_data中。
			{
				array_data[i][m * N_BYTE_INARRAY + j] = (record[m] & mask_groupby8[N_SLOTS - N_BYTE_INARRAY + j]) >> ((N_BYTE_INARRAY - 1 - j) * 8);
			}
		}
    }
	fclose(datafile);
}


/**
 * 将两个数的共同高位bit串，转换成SIMD数组，返回真实的共同bit个数.如果共同bit数大于length,返回length。
 * （注：这样做是因为，以后可能uint32_t的数只存储16位到SIMD，所以如果共同bit数不能大于N_BITS）
 */
void convert_bitcolumn_double(uint32_t left[N_ATTR], uint32_t right[N_ATTR], __m128i v[N_ATTR][N_BITS],int prefix_length[N_ATTR])
{
    int real_length = 0;
    int flag1, flag2;
    uint32_t set_value;
	for (int k = 0; k < N_ATTR; k++)
	{
		real_length = 0;
		for (int i = 0; i < N_BITS; i++)
		{
			flag1 = (left[k] >> (31 - i)) % 2;
			flag2 = (right[k] >> (31 - i)) % 2;
			if (flag1 != flag2)
			{
				break;

			}
			else {
				set_value = flag2 * UINT32MAX;
				v[k][i] = _mm_set_epi32(set_value, set_value, set_value, set_value);
				real_length++;
			}
		}
		prefix_length[k] = real_length;
	}
	//print2dmatrix(v);
}

/*
 tofind is the _m128i array which is converted by the value to find. The length
 indicates the actual length of SIMD array. 单点查询和范围查询均调用此函数
 */
void do_search(__m128i simd_mtx[N_ATTR][N_ROWS][N_BITS], __m128i to_find[N_ATTR][N_BITS], int prefix_length[N_ATTR])
{
    __m128i tmp;
	for (int k = 0; k < N_ATTR; k++)
	{
		tmp = _mm_setzero_si128();
		for (int row = 0; row < N_ROWS; row++)
		{
			for (int i = 0; i < prefix_length[k]; i++)
			{
				tmp = _mm_andnot_si128(_mm_xor_si128(simd_mtx[k][row][i], to_find[k][i]), ones); //只有当simd_mtx[row][i]和tofind[i]中对应的位相同时，tmp对应的位才是1.
				res_mask[row] = _mm_and_si128(res_mask[row], tmp);

				if (_mm_movemask_epi8(_mm_cmpeq_epi32(res_mask[row], zeros)) == 0xffff)
				{
					break;
				}
			}
		}
	}
	printf("\n");
	print1dmatrix(res_mask);
}


/**
 *检查数据是否在left和right之间，将位存储的数字组装成十进制的数字，然后和范围端点比较。
 real_length，将来用于，如果范围端点共同位少于DATATYPE_LEN-BITS_INARRAY时候的验证，暂未考虑这种情况
 */
void check(int prefix_length[N_ATTR], uint32_t left[N_ATTR], uint32_t right[N_ATTR])
{
    uint32_t pri_value[N_ATTR] ;   //记录要还原的数据
    uint32_t* res_mask_chunk;
    int is_one;
    uint32_t res_remain;
    int index = 0;
    int array_data_anchor_idx = 0;
	uint32_t newleft[N_ATTR];
	uint32_t newright[N_ATTR];
	for (int i = 0; i < N_ATTR; i++)
	{
		newleft[i] = left[i] & row_mask[N_BYTE_INARRAY];
		newright[i] = right[i] & row_mask[N_BYTE_INARRAY];
	}

   
	memset(pri_value, 0, sizeof(uint32_t)*N_ATTR);
    for (int i = 0; i < N_ROWS; i++)         //检查数组中的每个向量
    {
        if(_mm_movemask_epi8(_mm_cmpeq_epi32(res_mask[i], zeros)) == 0xffff)
        {
            continue;
        }
        
        res_mask_chunk = (uint32_t*)&res_mask[i];
        
        for (int j = 0; j < N_SLOTS; j++)    //检查向量中的每个槽
        {
            if (res_mask_chunk[j] == 0)
            {
                continue;
            }
            else
            {
                res_remain = res_mask_chunk[j];
                array_data_anchor_idx = (i << 7) + (j << 5) + 31;
                for (int k = 0; k < DATATYPE_LEN; k++)   //检查槽中的每个位
                {
					int flag = 1;
                    is_one = res_remain  & 0x00000001;  //从低位向高位check，得到数据的最低位
                    res_remain = res_remain >> 1;
                    if (is_one == 0)
                    {
                        continue;
                    }
                    else
                    {
                        index = array_data_anchor_idx - k; //通过i,j,k变量，定位出res_mask向量中为1位在array_data中的位置
						for (int m = 0; m < N_ATTR; m++)
						{
							int offset = m * N_BYTE_INARRAY;
							pri_value[m] = (array_data[index][0 + offset] << 16) + (array_data[index][1 + offset] << 8) + array_data[index][2 + offset];
							//printf("pri_value %d:%u ", m, pri_value[m]);
							if (pri_value[m] < newleft[m] || pri_value[m] > newright[m])
							{
								flag = 0;
								break;
							}
						}
						//printf("\n");
						if (flag == 1)
						{
							res_index++;
						}
                        //printf("pri_value:%d\t index=%d\n",pri_value,index);
                       
                        /**
                        else
                        {
                            res_mask_chunk[j] = res_mask_chunk[j] - power[k]; //将res_mask被验证无效的位置为0
                        }**/
                    }
                    /*if (res_remain == 0) //如果当前槽位剩下的数字都是0的话，停止移位，直接退出该槽位
                    {
                        break;
                    }*/
                }
            }
        }
    }
}



/**
* 范围查询
*/
void range_search(uint32_t left[N_ATTR], uint32_t right[N_ATTR])
{
	__m128i to_find[N_ATTR][N_BITS];
	int prefix_length[N_ATTR];
	convert_bitcolumn_double(left, right, to_find, prefix_length);
	do_search(matrix, to_find, prefix_length); //执行完,res_mask中即得到对应元素的掩码。
	check(prefix_length, left, right);
    
}


/**
 * 单点查询
 */
void single_search(uint32_t value[N_ATTR])
{
    range_search(value,value);
}


void clear_res_set()
{
    memset(res_set, 0, sizeof(uint32_t)*10000);
    res_index = 0;
}

int main()
{
    init(); //This cant be deleted
    
	laod_all_data(matrix);
	print3dmatrix(matrix);
	//print_array_data(array_data);

	//print2dmatrix(matrix);
	/*for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			printf("%u ", array_data[i][j]);
		}
		printf("\n");
	}*/
	uint32_t left[N_ATTR] = {4,4,4};
	uint32_t right[N_ATTR] = {9,9,9};
    clear_res_set();
	clock_t begin, end;
	begin = clock();
	//single_search(1);
    range_search(left, right);
    end = clock();
	printf("time with sse:%ld\n",(end - begin));
	printf("count:%d \n", res_index);
    
    /**
    for(int i=0; i<N; i++)
    {
        printf("index=%d\t%u\n",i,array_data[i][2]);
    }**/
    
    
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
