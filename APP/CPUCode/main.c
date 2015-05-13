#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <MaxSLiCInterface.h>
#include "Maxfiles.h"
#define DEBUG 1
#define EPS 1.0e-9
#define CUR_EPS 1.0e-9 

void max_print_matrix(double *m, int rows, int cols, char *hint_msg)
{
	fprintf(stdout, "====>%s\n", hint_msg);
	for(int i = 0; i < rows; i ++)
	{
		for(int j = 0; j < cols; j ++)
		{
			fprintf(stdout, "%g ", m[i * cols + j]);
		}
		fprintf(stdout, "\n");
	}
}
void jacobi_opt(double *A, double *x_base, double *b_base, int dim, int C, int total_equations, double *init_value, double *expected_error)
{
	int RUN = total_equations;
	double* x = malloc(dim*sizeof(double));
	double* z = malloc(dim*sizeof(double));
	double* b = b_base;
	int i, j, run;

	int stage = 0;
	for(run = 0; run < RUN; ++run) {
		double error = 1000;
		for(i = 0; i < dim; ++i) {
			x[i] = 0;
		}
		while(error > 1e-9) 
		{
			// update x
			for(i = 0; i < dim; ++i) {
				double sigma = 0;
				for(j = 0; j < dim; ++j) {
					if(j != i) {
						sigma += A[i*dim+j]*x[j];
					}
				}
				z[i] = (b[i]-sigma)/A[i*dim+i];
				++stage;
			}
			double *tmp = x;
			x = z;
			z = tmp;
			// compute error
			error = 0;
			for(i = 0; i < dim; ++i) {
				double y = 0;
				for(j = 0; j < dim; ++j) {
					y += A[i*dim+j]*x[j];
				}
				error += fabs(y-b[i]);
			}
		}
		expected_error[run] = error;
		memcpy(x_base+run*dim, x, dim*sizeof(double));

		// move to next problem to solve
		b += dim;
	}

	free(x);
	free(z);
}

void max_print_result(size_t dim_len, size_t equations, size_t MAX_ITER, double engine_time, double CPU_time )
{
	double speedup     = (fabs(engine_time) > 1e-9) ? (CPU_time/engine_time) : 0;
	const char *header = "Dim Len     Equations    MAX ITER    Engine Time    CPU Time    Speedup";
	char str_dim_len     [32] ;
	char str_equations   [32] ;
	char str_max_iter    [32] ;
	char str_engine_time [32] ;
	char str_CPU_time    [32] ;
	char str_speedup     [32] ;
	sprintf(str_dim_len     , "%ld"  , dim_len     );
	sprintf(str_equations   , "%ld"  , equations   );
	sprintf(str_max_iter    , "%ld"  , MAX_ITER    );
	sprintf(str_engine_time , "%.3f" , engine_time );
	sprintf(str_CPU_time    , "%.3f" , CPU_time    );
	sprintf(str_speedup     , "%.2f" , speedup     );

	char padding_space[32];
	memset(padding_space, 0, sizeof(padding_space));
	int  dim_len_padding = 0;
	dim_len_padding = strlen("Dim Len     ") - strlen(str_dim_len);
	for(int i = 0; i < dim_len_padding; i ++)
	{
			padding_space[i] = ' ';
	}
	strcat(str_dim_len, padding_space);
	
	memset(padding_space, 0, sizeof(padding_space));
	int equation_padding_len = 0;
	equation_padding_len = strlen("Equations    ") - strlen(str_equations);
	for(int i = 0; i < equation_padding_len; i ++)
	{
		padding_space[i] = ' ';
	}
	strcat(str_equations, padding_space);
	memset(padding_space, 0, sizeof(padding_space));
	int max_iter_padding_len = 0;
	max_iter_padding_len = strlen("MAX ITER    ") - strlen(str_max_iter);
	for(int i = 0; i < max_iter_padding_len; i ++)
	{
		padding_space[i] = ' ';
	}
	strcat(str_max_iter, padding_space);

	memset(padding_space, 0, sizeof(padding_space));
	int engine_time_padding_len = strlen("Engine Time    ") - strlen(str_engine_time);
	for(int i = 0; i < engine_time_padding_len; i ++)
	{
		padding_space[i] = ' ';
	}
	strcat(str_engine_time, padding_space);

	memset(padding_space, 0, sizeof(padding_space));
	int CPU_time_padding_len = strlen("CPU Time    ") - strlen(str_CPU_time);
	for(int i = 0; i < CPU_time_padding_len; i ++)
	{
		padding_space[i] = ' ';
	}
	strcat(str_CPU_time, padding_space);
	
	memset(padding_space, 0, sizeof(padding_space));
	int speedup_padding_len = strlen("CPU Time    ") - strlen(str_speedup);
	for(int i = 0; i < speedup_padding_len; i ++)
	{
		padding_space[i] = ' ';
	}
	strcat(str_speedup, padding_space);

	fprintf(stderr, "%s\n", header);
	fprintf(stderr, "%s%s%s%s%s%s\n",str_dim_len, str_equations,str_max_iter,str_engine_time,str_CPU_time,str_speedup);
}

void usage()
{
      printf("\nUsage:\n\n");
      printf("\t./JacobiRun [-d <number>] [-b <number >] [-i <number>]");
      printf("\nwhere:\n");
      printf("\t");
      printf("-d\tDemension Length, should between[2,200] and should be multipy of 2, default length is 2.\n");
      printf("\t");
      printf("-b\tBlocks to be caculated, each block contains 62 equations,should be\n\t\tan integer which is bigger than zero, default block number is 1.\n");
      printf("\t");
      printf("-i\tIteration times for Kernel, should bigger than 1, default iteration times is 20.\n");
      printf("\t");
      printf("-h\tPrint this help message.\n\n");
}

int main(int argc, char** argv) 
{
	max_file_t *max_file = jacobi_init();
	size_t dim               = 64;   // this should be a scalar input in the bitstream
	size_t MAX_ITER          = 20;
	size_t C                 = max_get_constant_uint64t(max_file, "C");
	size_t blks              = 100;
	size_t total_equations   = blks*C;
	clock_t engine_start     = 0;
	clock_t engine_end       = 0;
	double engine_total_time = 0.0;
	size_t max_dim           = max_get_constant_uint64t(max_file, "maxDimLen");
	if(argc == 1)
	{
		fprintf(stderr, "====>Info:Runing Jacobi with default parameter values:[Dimension = %ld, Iteration = %ld, blocks = %ld(%ld*%ld equations)], for details, see the README.txt\n", dim, MAX_ITER, blks, blks, C);
	}
	char *opt_str = "hd:b:i:";
	int opt = 0;
	int input_dim = dim;
	int input_iter = MAX_ITER;
	int input_blks = blks;
	while( (opt = getopt(argc, argv, opt_str)) != -1)
	{
		switch(opt)
		{
			case 'd':
				input_dim = atoi(optarg);	
				break;
			case 'b':
				input_blks = atoi(optarg);
				break;
			case 'i':
				input_iter = atoi(optarg);
				break;
			case 'h':
				usage();
				return 1;
			default:
				fprintf(stderr, "====>Error: Inputs contain invalid command line paramter(s)!\n");
				usage();
				return 1;
			}
		}
	max_file_free(max_file);
	if(input_dim <= 0 || input_dim > max_dim || input_dim % 2 != 0)
	{
		fprintf(stderr, "\n====>Error: Input dimension length is invalid, for details, see the usage below:\n");
		usage();
		return 1;
	}
	else
	{
		dim = (size_t)input_dim;
	}
	if(input_blks <= 0)
	{
		fprintf(stderr, "\n====>Error: Input block number is invalid, should bigger than zero.\n");
		usage();
		return 1;
	}
	else
	{
		blks = (size_t)input_blks;
	}
	if(input_iter <= 1)
	{
		fprintf(stderr, "\n====>Error: Input iteration number is invalid, should bigger than 1.\n");
		usage();
		return 1;
	}
	else
	{
		MAX_ITER = (size_t)input_iter;
	}
	total_equations = blks * C;

	double *A                  = malloc(dim*dim*sizeof(double));
	double *A_trans            = malloc(dim*dim*sizeof(double));
	double *b                  = malloc(total_equations*dim*sizeof(double));
	double *b_trans            = malloc(total_equations*dim*sizeof(double));
	double *diagA              = malloc(dim*sizeof(double));
	double *reverse_diagA      = malloc(dim*sizeof(double));
	double *x_init             = malloc(C*dim*sizeof(double));
	double *x_trans_init       = malloc(C*dim*sizeof(double));
	double *result             = malloc(total_equations * dim * sizeof(double));
	double *reorder_result     = malloc(total_equations * dim *sizeof(double));
	double *solutions          = malloc(total_equations * dim *sizeof(double));
	double *error              = malloc(total_equations*sizeof(double));
	double *error_bak          = malloc(total_equations*sizeof(double));
	int    *is_solution_valid  = malloc(total_equations*sizeof(int));
	int    *recacu_error_index = malloc(total_equations*sizeof(int));
	double *expected_error     = malloc(total_equations*sizeof(double));
	double *x_base             = malloc(total_equations * dim * sizeof(double));
	double *x_all_init         = malloc(total_equations * dim *sizeof(double));
	double *x_all_trans_init   = malloc(total_equations * dim *sizeof(double));
	memset(A,                0 , sizeof(double)*dim*dim);
	memset(A_trans,          0 , sizeof(double)*dim*dim);
	memset(b,                0 , sizeof(double)*dim*total_equations);
	memset(b_trans,          0 , sizeof(double)*dim*total_equations);
	memset(diagA,            0 , sizeof(double)*dim);
	memset(reverse_diagA,    0 , sizeof(double)*dim);
	memset(x_init,           0 , sizeof(double) *C*dim);
	memset(result,           0 , sizeof(double)*dim*total_equations);
	memset(reorder_result,   0 , sizeof(double)*dim*total_equations);
	memset(error,            0 , sizeof(double)*total_equations);
	memset(expected_error,   0 , sizeof(double)*total_equations);
	memset(x_base,           0 , sizeof(double)*dim*total_equations);
	memset(x_all_init,       0 , sizeof(double)*dim*total_equations);
	memset(x_all_trans_init, 0 , sizeof(double)*dim*total_equations);
	memset(is_solution_valid,0 , sizeof(int)*total_equations);
	
	for(int i = 0; i < total_equations; i ++)
	{
		recacu_error_index[i] = -1;
		expected_error[i]     = 1000;
		error_bak[i]          = 1000;
		for(int j = 0; j < dim; j ++)
		{
			solutions[i*dim + j] = 1000;
		}
	}

	/**
	 *  Generating random value for b and A
	 */
	srand(time(NULL));
	for(int i = 0; i < dim; ++i) {
		double sum = 0;
		for(int j = 0; j < dim; ++j) {
			if(i != j) {
				A[i*dim+j]     = 2.0*rand()/(double)RAND_MAX - 1 ; // random number between -1 and 1
				sum           += fabs(A[i*dim+j])                ;
			}
		}
		A[i * dim + i] = 1 + sum;
		diagA[i]       = 1.0/A[i * dim + i];
		reverse_diagA[i]  = A[i * dim + i];
	}
	
	double A_original[dim * dim];
	for(int i = 0; i < C*blks; i ++)
	{
			for(int j = 0; j < dim; j ++)
			{
				b[i * dim + j] = 2.0*rand()/(double)RAND_MAX - 1;
			}
	}

	for(int i = 0; i < dim; i ++)
	{
		for(int j = 0; j < dim; j ++)
		{
			A_original[i * dim + j] = A[i * dim + j];
			if(i != j)
			{
				A[i * dim + j] = A[i*dim + j] * diagA[i];
			}
		}
	}

	/**
	 * Reorder the input A and b 
	 */
	engine_start = clock();
	for(int i = 0; i < dim; i ++)
	{
		for(int j = 0; j < dim; j ++)
		{
			A_trans[i * dim + j] = A[j * dim + i];
		}
	}
	int count = 0;
	for(int yy = 0; yy < total_equations; yy += C)
	{
			for(int i = 0; i < dim; i ++)
			{
				for(int j = yy; j <yy + C; j ++)
				{
					b_trans[count] = b[j * dim + i]*diagA[i]; 
					count ++;
				}
			}
	}

	for(int k = 0; k < blks; k ++)
	{
			for ( int i = 0; i < C ; i ++ ) 
			{
					for ( int j = 0; j < dim; j ++ ) 
					{
							x_init[i * dim + j] = 0;
							x_trans_init[j*C + i] = x_init[i * dim + j];
					}
			}
		memcpy(x_all_trans_init + k * C * dim , x_trans_init , sizeof(double)*C*dim);
		memcpy(x_all_init       + k * C * dim , x_init       , sizeof(double)*C*dim);
	}
	
    jacobi(
		dim, 
		total_equations,
		MAX_ITER,
		A_trans                                ,
		dim * dim * sizeof(double)             ,
		b_trans                                ,
		total_equations * dim * sizeof(double) ,
		reverse_diagA                                  ,
		dim * sizeof(double)                   ,
		x_all_trans_init                       ,
		total_equations * dim * sizeof(double) ,
		error                                  ,
		total_equations * sizeof(double)       ,
		result                                 ,
	    total_equations * dim * sizeof(double) 
	  );

	for(int yy = 0; yy<total_equations; yy += C)
	{
			for(int i = 0; i <  C; i ++)
			{
				for(int j = 0; j < dim; j ++)
				{
					reorder_result[yy *dim + i*dim + j] = result[yy * dim + i + j * C];
				}
			}
	}

	/*Check Error to decide whether we need to restream into kernel again*/
	int recacu_cnt            = 0;
	int new_recacu_cnt        = 0;
	int actual_recacu_cnt     = 0;
	int new_actual_recacu_cnt = 0;

	double *x_latest_init       = malloc(total_equations * dim * sizeof(double)) ;
	double *x_latest_trans_init = malloc(total_equations * dim * sizeof(double)) ;
	double *recacu_b            = malloc(total_equations * dim * sizeof(double)) ;
	double *recacu_trans_b      = malloc(total_equations * dim * sizeof(double)) ;
	memset(x_latest_init       , 0 , total_equations * dim * sizeof(double))     ;
	memset(x_latest_trans_init , 0 , total_equations * dim * sizeof(double))     ;
	memset(recacu_b            , 0 , total_equations * dim * sizeof(double))     ;
	memset(recacu_trans_b      , 0 , total_equations * dim * sizeof(double))     ;

	int idx = 0;
	for(int i = 0; i < total_equations; i ++)
	{
		if(error[i] > CUR_EPS)
		{
			memcpy(x_latest_init + idx*dim, reorder_result + i*dim, dim*sizeof(double)); 
			memcpy(recacu_b      + idx*dim, b              + i*dim, dim*sizeof(double)); 
			recacu_error_index[idx] = i;			
			recacu_cnt ++        ;
			actual_recacu_cnt ++ ;
			idx ++;
		}
		else
		{
			error_bak[i] = error[i];
			memcpy(solutions +  i*dim, reorder_result + i*dim, dim*sizeof(double));
		}

	}
	while( recacu_cnt % C )
	{
		recacu_cnt ++;
	}

	/**
	 *  if recaculate count not zero, we start to restream data into kernel again
	 */
	int times = 1;
	while( recacu_cnt != 0 )
	{	
		/*Reorder Latest solutions init value */
		times ++;
		memset(x_latest_trans_init, 0, recacu_cnt*dim*sizeof(double));
		count = 0;
		for(int yy = 0; yy < recacu_cnt; yy += C)
		{
			for(int i = 0; i < dim; i ++)
			{
				for(int j = yy; j < yy + C; j ++)
				{
					x_latest_trans_init[count] = x_latest_init[j * dim + i]; 
					count ++;
				}
			}
		}


		/*Reorder latest b*/
		memset(recacu_trans_b, 0, total_equations*dim*sizeof(double));
		count = 0;
		for(int yy = 0; yy < recacu_cnt; yy += C)
		{
			for(int i = 0; i < dim; i ++)
			{
				for(int j = yy; j < yy + C; j ++)
				{
					recacu_trans_b[count] = recacu_b[j * dim + i]*diagA[i]; 
					count ++;
				}
			}
		}

		memset(error  , 0 , recacu_cnt * sizeof(double       )  ) ;
		memset(result , 0 , recacu_cnt * dim * sizeof(double )  ) ;
		jacobi(
			dim, 
			recacu_cnt,
			MAX_ITER,
			A_trans                           ,
			dim * dim * sizeof(double)        ,
			recacu_trans_b                    ,
			recacu_cnt * dim * sizeof(double) ,
			reverse_diagA                             ,
			dim * sizeof(double)              ,
			x_latest_trans_init               ,
			recacu_cnt * dim * sizeof(double) ,
			error                             ,
			recacu_cnt * sizeof(double)       ,
			result                            ,
			recacu_cnt * dim * sizeof(double)
		  );


		for(int yy = 0; yy < recacu_cnt; yy += C)
		{
			for(int i = 0; i <  C; i ++)
			{
				for(int j = 0; j < dim; j ++)
				{
					reorder_result[yy *dim + i*dim + j] = result[yy * dim + i + j * C];
				}
			}
		}
		
		new_recacu_cnt = 0;
		new_actual_recacu_cnt = 0;
		int idx2 = 0;
		for(int i = 0; i < actual_recacu_cnt; i ++)
		{
			if(error[i] > CUR_EPS)
			{
				memcpy(x_latest_init + new_recacu_cnt*dim, reorder_result + i*dim, dim*sizeof(double)); 
				memcpy(recacu_b      + new_recacu_cnt*dim, recacu_b       + i*dim, dim*sizeof(double)); 
				recacu_error_index[idx2] = recacu_error_index[i]; 
				new_recacu_cnt ++;
				new_actual_recacu_cnt ++;
				idx2 ++;
			}
			else
			{
				error_bak[ recacu_error_index[i]] = error[i];
				memcpy(solutions + recacu_error_index[i] *dim, reorder_result + i*dim, dim*sizeof(double));
			}
		}

		/* padding to multipy of C */
		while( new_recacu_cnt % C )
		{
			new_recacu_cnt ++;
		}
		/* update the current recaculating solution numbers */
		recacu_cnt = new_recacu_cnt;
		actual_recacu_cnt = new_actual_recacu_cnt;
	}//loop while

	engine_end        = clock();
	engine_total_time = (double)(engine_end - engine_start) / CLOCKS_PER_SEC;
	fprintf(stderr, "=========>Kernel Complete, Stream Times: %d\n", times);
	clock_t cpu_start = clock();
	jacobi_opt(A_original, x_base, b, dim, C, total_equations, x_all_init , expected_error);
	clock_t cpu_end = clock();
	double cpu_total_time = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

	/* Compare the result with the standard result */
	int cnt = 0;
	int index = 0;
	for(int i = 0; i < total_equations; i ++)
	{
		for(int j = 0; j < dim; j ++)
		{
			double diff = solutions[i * dim + j] - x_base[i*dim + j];
			if(fabs(diff) > EPS)
			{
					fprintf(stderr, "error: atual=%.10f, expect=%.10f, err=%.10e\n",
							solutions[i * dim + j], x_base[i*dim + j], diff);
					cnt ++;
					index ++;
			}
		}
	}
	if(cnt == 0)
	{
		max_print_result(dim, total_equations, MAX_ITER, engine_total_time, cpu_total_time);
		fprintf(stderr, "==========>All Test Passed\n\n");
	}
	else
	{
		fprintf(stderr, "!!!Test Failed:%d\n\n", cnt);
	}

	free ( A                   ) ;
	free ( A_trans             ) ;
	free ( b                   ) ;
	free ( b_trans             ) ;
	free ( diagA               ) ;
	free ( reverse_diagA       ) ;
	free ( x_init              ) ;
	free ( error               ) ;
	free ( error_bak           ) ;
	free ( recacu_error_index  ) ;
	free ( expected_error      ) ;
	free ( result              ) ;
	free ( reorder_result      ) ;
	free ( solutions           ) ;
	free ( x_base              ) ;
	free ( x_all_init          ) ;
	free ( x_all_trans_init    ) ;
	free ( x_latest_init       ) ;
	free ( x_latest_trans_init ) ;
	free ( recacu_b            ) ;

	int status = (cnt == 0) ? 0:1;
	return status;
}
