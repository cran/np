/* Copyright (C) J. Racine, 1995-2001 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <assert.h>

#include "headers.h"
#include "matrix.h"

#ifdef MPI

#include "mpi.h"

extern  int my_rank;
extern  int source;
extern  int dest;
extern  int tag;
extern  int iNum_Processors;
extern  int iSeed_my_rank;
extern  MPI_Status status;
#endif

/*
int int_LARGE_SF; 
int int_DEBUG;
int int_VERBOSE;
int int_NOKEYPRESS;
int int_DISPLAY_CV;
int int_RANDOM_SEED;
int int_MINIMIZE_IO;
int int_ORDERED_CATEGORICAL_GRADIENT;
int int_PREDICT; 
int int_ROBUST;
int int_SIMULATION;
int int_TAYLOR;
int int_WEIGHTS;
*/

#ifdef RCSID
static char rcsid[] = "$Id: jksum.c,v 1.16 2006/11/02 16:56:49 tristen Exp $";
#endif

/* Some externals for numerical routines */

extern int num_obs_train_extern;
extern int num_obs_eval_extern;
extern int num_var_continuous_extern;
extern int num_var_unordered_extern;
extern int num_var_ordered_extern;
extern int num_reg_continuous_extern;
extern int num_reg_unordered_extern;
extern int num_reg_ordered_extern;
extern int *num_categories_extern;
extern double **matrix_categorical_vals_extern;

extern double **matrix_X_continuous_train_extern;
extern double **matrix_X_unordered_train_extern;
extern double **matrix_X_ordered_train_extern;
extern double **matrix_X_continuous_eval_extern;
extern double **matrix_X_unordered_eval_extern;
extern double **matrix_X_ordered_eval_extern;

extern double **matrix_Y_continuous_train_extern;
extern double **matrix_Y_unordered_train_extern;
extern double **matrix_Y_ordered_train_extern;
extern double **matrix_Y_continuous_eval_extern;
extern double **matrix_Y_unordered_eval_extern;
extern double **matrix_Y_ordered_eval_extern;

extern double *vector_Y_extern;
extern double *vector_T_extern;
extern double *vector_Y_eval_extern;

/* Quantile - no Y ordered or unordered used, but defined anyways */

extern double **matrix_Y_continuous_quantile_extern;
extern double **matrix_Y_unordered_quantile_extern;
extern double **matrix_Y_ordered_quantile_extern;
extern double **matrix_X_continuous_quantile_extern;
extern double **matrix_X_unordered_quantile_extern;
extern double **matrix_X_ordered_quantile_extern;

extern int int_ll_extern;

extern int KERNEL_reg_extern;
extern int KERNEL_reg_unordered_extern;
extern int KERNEL_reg_ordered_extern;
extern int KERNEL_den_extern;
extern int KERNEL_den_unordered_extern;
extern int KERNEL_den_ordered_extern;
extern int BANDWIDTH_reg_extern;
extern int BANDWIDTH_den_extern;

extern int itmax_extern;
extern double small_extern;

/* Statics for dependence metric */

extern int num_lag_extern;
extern int int_lag_extern;
extern int int_iter_extern;

extern double *vector_scale_factor_dep_met_bivar_extern;
extern double *vector_scale_factor_dep_met_univar_extern;
extern double *vector_scale_factor_dep_met_univar_lag_extern;

extern double y_min_extern;
extern double y_max_extern;

int kernel_convolution_weighted_sum(
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_reg,
int num_obs_train,
int num_obs_eval,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double *vector_Y,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *kernel_sum)
{

	/* This function takes a vector Y and returns a convolution kernel
		  weighted sum. By default Y should be a vector of ones (simply
		  compute the kernel sum). This function will allow users to `roll
		  their own' with mixed data convolution kernel sums. */

	/* Declarations */

	int i;
	int j;
	int l;

	double prod_kernel;
	double sum_y_ker;

	double *lambda;
	double **matrix_bandwidth = NULL;

	double *psum;
	double *py;

#ifdef MPI
	int stride = ceil((double) num_obs_eval / (double) iNum_Processors);
	if(stride < 1) stride = 1;
#endif

	/* Allocate memory for objects */

	lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);

	if((BANDWIDTH_reg == 0)||(BANDWIDTH_reg == 1))
	{
		matrix_bandwidth = alloc_matd(num_obs_eval,num_reg_continuous);
	}
	else if(BANDWIDTH_reg == 2)
	{
		matrix_bandwidth = alloc_matd(num_obs_train,num_reg_continuous);
	}

	/* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

	if(kernel_bandwidth_mean(
		KERNEL_reg,
		BANDWIDTH_reg,
		num_obs_train,
		num_obs_eval,
		0,
		0,
		0,
		num_reg_continuous,
		num_reg_unordered,
		num_reg_ordered,
		vector_scale_factor,
		matrix_X_continuous_train,	 /* Not used */
		matrix_X_continuous_eval,		 /* Not used */
		matrix_X_continuous_train,
		matrix_X_continuous_eval,
		matrix_bandwidth,						 /* Not used */
		matrix_bandwidth,
		lambda) == 1)
	{
#ifndef MPI
		printf("\n** Error: invalid bandwidth.");
		printf("\nProgram Terminated.\n");
		exit(EXIT_FAILURE);
#endif
#ifdef MPI
		if(my_rank == 0)
		{
			printf("\n** Error: invalid bandwidth.");
			printf("\nProgram Terminated.\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(EXIT_FAILURE);
#endif
	}

#ifndef MPI

	if(BANDWIDTH_reg == 0)
	{

		psum = &kernel_sum[0];

		for(j=0; j < num_obs_eval; j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0],
						matrix_bandwidth[l][0],
						matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			*psum++ = sum_y_ker;

		}

	}
	else if(BANDWIDTH_reg == 1)
	{

		psum = &kernel_sum[0];

		for(j=0; j < num_obs_eval; j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j],
						matrix_bandwidth[l][i],
						matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			*psum++ = sum_y_ker;

		}

	}
	else
	{

		psum = &kernel_sum[0];

		for(j=0; j < num_obs_eval; j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i],
						matrix_bandwidth[l][j],
						matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			*psum++ = sum_y_ker;

		}

	}
#endif

#ifdef MPI

	if(BANDWIDTH_reg == 0)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][0],
						matrix_bandwidth[l][0],
						matrix_bandwidth[l][0]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			kernel_sum[j-my_rank*stride] = sum_y_ker;

		}

	}
	else if(BANDWIDTH_reg == 1)
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][j],
						matrix_bandwidth[l][i],
						matrix_bandwidth[l][j]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			kernel_sum[j-my_rank*stride] = sum_y_ker;

		}

	}
	else
	{

		for(j=my_rank*stride; (j < num_obs_eval) && (j < (my_rank+1)*stride); j++)
		{

			sum_y_ker = 0.0;
			py = &vector_Y[0];

			for(i=0; i < num_obs_train; i++)
			{

				prod_kernel = 1.0;

				for(l = 0; l < num_reg_continuous; l++)
				{
					prod_kernel *= kernel_convol(KERNEL_reg,BANDWIDTH_reg,
						(matrix_X_continuous_eval[l][j]-matrix_X_continuous_train[l][i])/matrix_bandwidth[l][i],
						matrix_bandwidth[l][j],
						matrix_bandwidth[l][i]);
				}

				for(l = 0; l < num_reg_unordered; l++)
				{
					prod_kernel *= kernel_unordered_convolution(KERNEL_unordered_reg,
						matrix_X_unordered_eval[l][j],
						matrix_X_unordered_train[l][i],
						lambda[l],
						num_categories[l],
						matrix_categorical_vals[l]);
				}

				for(l = 0; l < num_reg_ordered; l++)
				{
					prod_kernel *= kernel_ordered_convolution(KERNEL_ordered_reg,
						matrix_X_ordered_eval[l][j],
						matrix_X_ordered_train[l][i],
						lambda[l+num_reg_unordered],
						num_categories[l+num_reg_unordered],
						matrix_categorical_vals[l+num_reg_unordered]);
				}

				sum_y_ker +=  *py++ *prod_kernel;

			}

			kernel_sum[j-my_rank*stride] = sum_y_ker;

		}

	}

	MPI_Gather(kernel_sum, stride, MPI_DOUBLE, kernel_sum, stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(kernel_sum, num_obs_eval, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	free(lambda);

	free_mat(matrix_bandwidth,num_reg_continuous);

	return(0);

}

#define ONE_OVER_SQRT_TWO_PI 0.39894228040143267794

double np_gauss2(const double z){
  return ONE_OVER_SQRT_TWO_PI*exp(-0.5*z*z);
}

double np_gauss4(const double z){
  return ONE_OVER_SQRT_TWO_PI*(1.5-0.5*z*z)*exp(-0.5*z*z);
}

double np_gauss6(const double z){
  const double z2 = z*z;
  return ONE_OVER_SQRT_TWO_PI*(1.875-z2*(1.25+z2*0.125))*exp(-0.5*z2);
}

double np_gauss8(const double z){
  const double z2 = z*z;
  return ONE_OVER_SQRT_TWO_PI*(2.1875-z2*(2.1875+z2*(0.4375-z2*0.02083333333)))*exp(-0.5*z2);
}

double np_epan2(const double z){
  return (z*z < 5.0)?((double)(0.33541019662496845446-0.067082039324993690892*z*z)):0.0;
}

double np_epan4(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?((double)(0.008385254916*(-15.0+7.0*z2)*(-5.0+z2))):0.0;
}

double np_epan6(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?((double)(0.33541019662496845446*(2.734375-3.28125*z2+0.721875*z2)*(1.0-0.2*z2))):0.0;
}

double np_epan8(const double z){
  const double z2 = z*z;
  return (z2 < 5.0)?((double)(0.33541019662496845446*(3.5888671875-7.8955078125*z2
                                                      +4.1056640625*z2-.5865234375*z2)
                              *(1.0-0.2*z2))):0.0;
}

double np_rect(const double z){
  return (z*z < 1.0)?0.5:0.0;
}

double np_uaa(const int same_cat,const double lambda, const int c){
  return (same_cat)?(1.0-lambda):lambda/((double)c-1.0);
}

double np_uli_racine(const int same_cat, const double lambda, const int c){
  return (same_cat)?1.0:lambda;
}

double np_owang_van_ryzin(const double x, const double y, const double lambda){
  return (x == y)?(1.0-lambda):ipow(lambda, (int)fabs(x-y))*(1.0-lambda)*0.5;
}

double np_oli_racine(const double x, const double y, const double lambda){
  return (x == y)?1.0:ipow(lambda, (int)fabs(x-y));

}

/* convolution kernels */

/* should be added */


double (* const allck[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, 
                                  np_epan2, np_epan4, np_epan6, np_epan8, 
                                  np_rect };
double (* const allok[])(double, double, double) = { np_owang_van_ryzin, np_oli_racine };
double (* const alluk[])(int, double, int) = { np_uaa, np_uli_racine };


/* 
   np_kernelv does weighted products of vectors - this is useful for 
   product kernels, where each kernel in each dimension acts as a weight.
*/

/* xt = training data */
/* xw = x weights */

void np_ckernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double h, 
                 double * const result){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i,j; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, 
                                   np_epan2, np_epan4, np_epan6, np_epan8, 
                                   np_rect };

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*k[KERNEL]((xt[i]-x)/h);

}

void np_ukernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double lambda, const int ncat,
                 double * const result){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(int, double, int) = { np_uaa, np_uli_racine };

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*k[KERNEL]((xt[i]==x), lambda, ncat);

}

void np_okernelv(const int KERNEL, 
                 const double * const xt, const int num_xt, 
                 const int do_xw,
                 const double x, const double lambda,
                 double * const result){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i; 
  int j, bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double (* const k[])(double, double, double) = { np_owang_van_ryzin, np_oli_racine };

  for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw)
    result[i] = xw[j]*k[KERNEL](xt[i], x, lambda);
}


void np_outer_weighted_sum(double * const * const mat_A, const int ncol_A, 
                           double * const * const mat_B, const int ncol_B,
                           double * const weights, const int num_weights,
                           const int do_leave_one_out, const int leave_one_out,
                           const int kpow,
                           const int parallel_sum,
                           const int bandwidth_divide, const double dband,
                           double * const result){

  int i,j,k;
  const int kstride = (parallel_sum ? (MAX(ncol_A, 1)*MAX(ncol_B, 1)) : 0);
  const int max_A = MAX(ncol_A, 1);
  const int max_B = MAX(ncol_B, 1);
  const int have_A = (ncol_A == 0 ? 0 : 1);
  const int have_B =  (ncol_B == 0 ? 0 : 1);

  double unit_weight = 1.0;
  double * const punit_weight = &unit_weight;

  double * const * const pmat_A = (ncol_A == 0 ? 0 : 1)?
    mat_A:&punit_weight;

  double * const * const pmat_B = (ncol_B == 0 ? 0 : 1)?
    mat_B:&punit_weight;

  const double db = (bandwidth_divide ? dband : unit_weight);
  double temp = FLT_MAX;

  if (do_leave_one_out) {
    temp = weights[leave_one_out];
    weights[leave_one_out] = 0.0;
  }

  for (k = 0; k < num_weights; k++)
    for (j = 0; j < max_A; j++)
      for (i = 0; i < max_B; i++)
        result[k*kstride+j*max_B+i] += pmat_A[j][k*have_A]*pmat_B[i][k*have_B]*ipow(weights[k]/db, kpow);

  if (do_leave_one_out)
    weights[leave_one_out] = temp;
}

  /* 
     np_kernelv_sum does weighted product sums of vectors - this is useful for 
     product kernels, where each kernel in each dimension acts as a weight.
  */
  /*
    xw1 is applied before exponentiation, xw2 is applied subsequently
  */

double np_ckernelv_sum_double_weighted(const int KERNEL, 
                                       double * xt, int num_xt, 
                                       double * xw1, double *xw2, 
                                       double x, double h, int kpow
                                       ){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i;
  double ksum = 0.0;

  double (* const k[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, 
                                   np_epan2, np_epan4, np_epan6, np_epan8, 
                                   np_rect };
  if(xw1 != NULL)
    if(xw2 != NULL)
      for (i = 0; i < num_xt; i++)
        ksum += ipow(xw1[i]*k[KERNEL]((xt[i]-x)/h),kpow)*xw2[i];
    else
      for (i = 0; i < num_xt; i++)
        ksum += ipow(xw1[i]*k[KERNEL]((xt[i]-x)/h),kpow);
  else
    if(xw2 != NULL)
      for (i = 0; i < num_xt; i++)
        ksum += ipow(k[KERNEL]((xt[i]-x)/h),kpow)*xw2[i];
    else
      for (i = 0; i < num_xt; i++)
        ksum += ipow(k[KERNEL]((xt[i]-x)/h),kpow);


  return ksum;
}

void np_ckernelv_psum_double_weighted(const int KERNEL, 
                                      double * xt, int num_xt, 
                                      double * xw1, double xw2, 
                                      double x, double h, int kpow,
                                      double * result
                                      ){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i;

  double (* const k[])(double) = { np_gauss2, np_gauss4, np_gauss6, np_gauss8, 
                                   np_epan2, np_epan4, np_epan6, np_epan8, 
                                   np_rect };
  if(xw1 != NULL)
    for (i = 0; i < num_xt; i++)
      result[i] += ipow(xw1[i]*k[KERNEL]((xt[i]-x)/h),kpow)*xw2;
  else
    for (i = 0; i < num_xt; i++)
      result[i] += ipow(k[KERNEL]((xt[i]-x)/h),kpow)*xw2;

}

double np_ukernelv_sum_double_weighted(const int KERNEL, 
                                       double * xt, int num_xt, 
                                       double * xw1, double *xw2, 
                                       double x, double lambda, int ncat,
                                       int kpow
                                       ){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i;
  double ksum = 0.0;

  double (* const k[])(int, double, int) = { np_uaa, np_uli_racine };

  if(xw1 != NULL)
    if(xw2 != NULL)
      for (i = 0; i < num_xt; i++)
        ksum += ipow(xw1[i]*k[KERNEL]((xt[i]==x), lambda, ncat), kpow)*xw2[i];
    else
      for (i = 0; i < num_xt; i++)
        ksum += ipow(xw1[i]*k[KERNEL]((xt[i]==x), lambda, ncat), kpow);
  else
    if(xw2 != NULL)
      for (i = 0; i < num_xt; i++)
        ksum += ipow(k[KERNEL]((xt[i]==x), lambda, ncat), kpow)*xw2[i];
    else
      for (i = 0; i < num_xt; i++)
        ksum += ipow(k[KERNEL]((xt[i]==x), lambda, ncat), kpow);

  return ksum;
}

void np_ukernelv_psum_double_weighted(const int KERNEL, 
                                      double * xt, int num_xt, 
                                      double * xw1, double xw2, 
                                      double x, double lambda, int ncat,
                                      int kpow, double * result
                                      ){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i;

  double (* const k[])(int, double, int) = { np_uaa, np_uli_racine };

  if(xw1 != NULL)
    for (i = 0; i < num_xt; i++)
      result[i] += ipow(xw1[i]*k[KERNEL]((xt[i]==x), lambda, ncat), kpow)*xw2;
  else
    for (i = 0; i < num_xt; i++)
      result[i] += ipow(k[KERNEL]((xt[i]==x), lambda, ncat), kpow)*xw2;

}

double np_okernelv_sum_double_weighted(const int KERNEL, 
                                       double * xt, int num_xt, 
                                       double * xw1, double *xw2, 
                                       double x, double lambda,
                                       int kpow
                                       ){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i;
  double ksum = 0.0;

  double (* const k[])(double, double, double) = { np_owang_van_ryzin, np_oli_racine };

  if(xw1 != NULL)
    if(xw2 != NULL)
      for (i = 0; i < num_xt; i++)
        ksum += ipow(xw1[i]*k[KERNEL](xt[i], x, lambda), kpow)*xw2[i];
    else
      for (i = 0; i < num_xt; i++)
        ksum += ipow(xw1[i]*k[KERNEL](xt[i], x, lambda), kpow);
  else
    if(xw2 != NULL)
      for (i = 0; i < num_xt; i++)
        ksum += ipow(k[KERNEL](xt[i], x, lambda), kpow)*xw2[i];
    else
      for (i = 0; i < num_xt; i++)
        ksum += ipow(k[KERNEL](xt[i], x, lambda), kpow);

  return ksum;
}

void np_okernelv_psum_double_weighted(const int KERNEL, 
                                      double * xt, int num_xt, 
                                      double * xw1, double xw2, 
                                      double x, double lambda,
                                      int kpow, double *result
                                      ){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i;

  double (* const k[])(double, double, double) = { np_owang_van_ryzin, np_oli_racine };

  if(xw1 != NULL)
    for (i = 0; i < num_xt; i++)
      result[i] += ipow(xw1[i]*k[KERNEL](xt[i], x, lambda), kpow)*xw2;
  else
    for (i = 0; i < num_xt; i++)
      result[i] += ipow(k[KERNEL](xt[i], x, lambda), kpow)*xw2;

}


#define BW_FIXED   0
#define BW_GEN_NN  1
#define BW_ADAP_NN 2

int kernel_weighted_sum_np(
const int KERNEL_reg,
const int KERNEL_unordered_reg,
const int KERNEL_ordered_reg,
const int BANDWIDTH_reg,
const int num_obs_train,
const int num_obs_eval,
const int num_reg_unordered,
const int num_reg_ordered,
const int num_reg_continuous,
const int leave_one_out,
const int kernel_pow,
const int bandwidth_divide,
int do_smooth_coef_weights,
double * const * const matrix_X_unordered_train,
double **matrix_X_ordered_train,
double **matrix_X_continuous_train,
double **matrix_X_unordered_eval,
double **matrix_X_ordered_eval,
double **matrix_X_continuous_eval,
double **matrix_Y,
double **matrix_W,
double *vector_scale_factor,
int *num_categories,
double *weighted_sum){
  
  /* This function takes a vector Y and returns a kernel weighted
     leave-one-out sum. By default Y should be a vector of ones
     (simply compute the kernel sum). This function will allow users
     to `roll their own' with mixed data leave-one out kernel sums. */

  /* num_var_ordered_extern contains number of columns of the weight matrix */


  /* Declarations */

  int i,j,l, mstep, js, je, num_obs_eval_alloc, sum_element_length;
  int do_psum; 

  const int num_xt = (BANDWIDTH_reg == BW_ADAP_NN)?num_obs_eval:num_obs_train;
  const int ws_step = (BANDWIDTH_reg == BW_ADAP_NN)? 0 :
    (MAX(num_var_continuous_extern, 1) * MAX(num_var_ordered_extern, 1));

  double *lambda, **matrix_bandwidth, *m = NULL;
  double *tprod, dband, *buf, *ws;

  double * const * const xtc = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_continuous_eval:matrix_X_continuous_train;
  double * const * const xtu = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_unordered_eval:matrix_X_unordered_train;
  double * const * const xto = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_ordered_eval:matrix_X_ordered_train;

  double * const * const xc = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_continuous_train:matrix_X_continuous_eval;
  double * const * const xu = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_unordered_train:matrix_X_unordered_eval;
  double * const * const xo = (BANDWIDTH_reg == BW_ADAP_NN)?
    matrix_X_ordered_train:matrix_X_ordered_eval;

#ifdef MPI
  int stride = MAX(ceil((double) num_obs_eval / (double) iNum_Processors),1);
  num_obs_eval_alloc = stride*iNum_Processors;
#else
  num_obs_eval_alloc = num_obs_eval;
#endif

  if (num_obs_eval == 0) {
    return(1);
  }

  do_psum = BANDWIDTH_reg == BW_ADAP_NN;
  /* Allocate memory for objects */

  mstep = (BANDWIDTH_reg==BW_GEN_NN)?num_obs_eval:
    ((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_train:1);

  lambda = alloc_vecd(num_reg_unordered+num_reg_ordered);
  matrix_bandwidth = alloc_tmatd(mstep, num_reg_continuous);  
  
  tprod = alloc_vecd((BANDWIDTH_reg==BW_ADAP_NN)?num_obs_eval:num_obs_train);

  sum_element_length = MAX(num_var_continuous_extern, 1) * 
    MAX(num_var_ordered_extern, 1);

  buf = alloc_vecd(num_obs_eval_alloc*sum_element_length);

  /* clear memory in temporary buffer and in result */
  for(i = 0; i < num_obs_eval_alloc*sum_element_length; i++){
    buf[i] = 0.0;
    weighted_sum[i] = 0.0;
  }

  /* assert(!(BANDWIDTH_reg == BW_ADAP_NN)); */
  /* Conduct the estimation */

  /* Generate bandwidth vector given scale factors, nearest neighbors, or lambda */

  if(kernel_bandwidth_mean(
                           KERNEL_reg,
                           BANDWIDTH_reg,
                           num_obs_train,
                           num_obs_eval,
                           0,
                           0,
                           0,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           vector_scale_factor,
                           matrix_X_continuous_train,				 /* Not used */
                           matrix_X_continuous_eval,				 /* Not used */
                           matrix_X_continuous_train,
                           matrix_X_continuous_eval,
                           matrix_bandwidth,						 /* Not used */
                           matrix_bandwidth,
                           lambda)==1){

    free(lambda);
    free_tmat(matrix_bandwidth);
    free(tprod);

    return(1);
  }

  if ((num_obs_train != num_obs_eval) && leave_one_out){
    
    printf("\ntraining and evaluation data must be the same to use leave one out estimator");
    free(lambda);
    free_tmat(matrix_bandwidth);
    free(tprod);
    return(1);
  }

  if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
#ifdef MPI
    js = stride * my_rank;
    je = MIN(num_obs_eval - 1, js + stride - 1);
#else
    js = 0;
    je = num_obs_eval - 1;
#endif
    
    ws = weighted_sum + js * sum_element_length;
  } else {
#ifdef MPI
    js = stride * my_rank;
    je = MIN(num_obs_train - 1, js + stride - 1);
    ws = buf;
#else
    js = 0;
    je = num_obs_train - 1;
    ws = weighted_sum;
#endif
  }
  
    /* do sums */
  for(j=js; j <= je; j++, ws += ws_step){
    dband = 1.0;

    if (num_reg_continuous > 0){
      m = matrix_bandwidth[0];
      if (BANDWIDTH_reg != BW_FIXED)
        m += j;
    }

    l = 0;

    /* continuous first */

    /* for the first iteration, no weights */
    /* for the rest, the accumulated products are the weights */
    
    for(i=0; i < num_reg_continuous; i++, l++, m += mstep){
      np_ckernelv(KERNEL_reg, xtc[i], num_xt, l, xc[i][j], *m, tprod);
      dband *= *m;
    }

    /* unordered second */

    for(i=0; i < num_reg_unordered; i++, l++){
      np_ukernelv(KERNEL_unordered_reg, xtu[i], num_xt, l, xu[i][j], 
                  lambda[i], num_categories[i], tprod);
    }

    /* ordered third */
    for(i=0; i < num_reg_ordered; i++, l++){
      np_okernelv(KERNEL_ordered_reg, xto[i], num_xt, l,
                  xo[i][j], lambda[num_reg_unordered+i], 
                  tprod);
    }

    /* expand matrix outer product, multiply by kernel weights, etc, do sum */

    np_outer_weighted_sum(matrix_W, num_var_ordered_extern, 
                          matrix_Y, num_var_continuous_extern,
                          tprod, num_obs_train,
                          leave_one_out, j,
                          kernel_pow,
                          do_psum,
                          bandwidth_divide, dband,
                          ws);
  }

#ifdef MPI
  if (BANDWIDTH_reg == BW_FIXED || BANDWIDTH_reg == BW_GEN_NN){
    MPI_Allgather(weighted_sum + js * sum_element_length, stride * sum_element_length, MPI_DOUBLE, weighted_sum, stride * sum_element_length, MPI_DOUBLE, MPI_COMM_WORLD);
  } else if(BANDWIDTH_reg == BW_ADAP_NN){
    MPI_Allreduce(ws, weighted_sum, num_obs_eval*sum_element_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  free(lambda);
  free_tmat(matrix_bandwidth);
  free(tprod);
  free(buf);
  
  return(0);
}


int np_kernel_estimate_con_density_categorical_convolution_cv(
int KERNEL_den,
int KERNEL_unordered_den,
int KERNEL_ordered_den,
int KERNEL_reg,
int KERNEL_unordered_reg,
int KERNEL_ordered_reg,
int BANDWIDTH_den,
int num_obs,
int num_var_unordered,
int num_var_ordered,
int num_var_continuous,
int num_reg_unordered,
int num_reg_ordered,
int num_reg_continuous,
double **matrix_Y_unordered,
double **matrix_Y_ordered,
double **matrix_Y_continuous,
double **matrix_X_unordered,
double **matrix_X_ordered,
double **matrix_X_continuous,
double *vector_scale_factor,
int *num_categories,
double **matrix_categorical_vals,
double *cv){


  int i, j, k, l, m, n, o, ib, ob, ms, ns, os;

  const int blklen = MIN(64, num_obs);
  int blj, bli, blk;

  double * sum_ker_convol, * sum_ker_marginal, * sum_ker;

  double *lambda;
  double **matrix_bandwidth_var, **matrix_bandwidth_reg;

  double * blk_xi, * blk_xj, * blk_yij, ts;

  double (* const yck)(double) = allck[KERNEL_den];
  double (* const xck)(double) = allck[KERNEL_reg];
  double (* const yok)(double, double, double) = allok[KERNEL_ordered_den];
  double (* const xok)(double, double, double) = allok[KERNEL_ordered_reg];
  double (* const yuk)(int, double, int) = alluk[KERNEL_unordered_den];
  double (* const xuk)(int, double, int) = alluk[KERNEL_unordered_reg];


  lambda = alloc_vecd(num_var_unordered+num_reg_unordered+num_var_ordered+num_reg_ordered);
  matrix_bandwidth_var = alloc_matd(num_obs,num_var_continuous);
  matrix_bandwidth_reg = alloc_matd(num_obs,num_reg_continuous);

  blk_xi = (double *)malloc(sizeof(double)*blklen*blklen);
  assert(blk_xi != NULL);

  blk_xj = (double *)malloc(sizeof(double)*blklen*blklen);
  assert(blk_xj != NULL);

  blk_yij = (double *)malloc(sizeof(double)*blklen*blklen);
  assert(blk_xj != NULL);

  sum_ker = (double *)calloc(num_obs, sizeof(double));
  assert(sum_ker != NULL);

  sum_ker_convol = (double *)calloc(num_obs, sizeof(double));
  assert(sum_ker_convol != NULL);

  sum_ker_marginal = (double *)calloc(num_obs, sizeof(double));
  assert(sum_ker_marginal != NULL);

  assert(BANDWIDTH_den == BW_FIXED);

  if(kernel_bandwidth_mean(KERNEL_den,
                           BANDWIDTH_den,
                           num_obs,
                           num_obs,
                           num_var_continuous,
                           num_var_unordered,
                           num_var_ordered,
                           num_reg_continuous,
                           num_reg_unordered,
                           num_reg_ordered,
                           vector_scale_factor,
                           matrix_Y_continuous,
                           matrix_Y_continuous,
                           matrix_X_continuous,
                           matrix_X_continuous,
                           matrix_bandwidth_var,
                           matrix_bandwidth_reg,
                           lambda)==1){
    free(lambda);
    free_mat(matrix_bandwidth_var,num_var_continuous);
    free_mat(matrix_bandwidth_reg,num_reg_continuous);
    return(1);
  }

  *cv = 0.0;

  // top level loop corresponds to a block of X_l's 
  // (leave one out evaluation points)
  for(k=0; k < num_obs; k+=blklen){ 
    blk = MIN(num_obs-blklen, k);

    // we generate blocks of X_j, X_i, and Y_ji
    // outermost loop determines the X_j's
    // innermost loop we generate new X_i and Y_ji blocks

    // one improvement would be a flip-flop algorithm where the roles
    // of X_j and X_i are swapped and only a new Y_ji need be generated
    for(j=0; j < num_obs; j+=blklen){
      blj = MIN(num_obs-blklen, j);

      // first: fill in blk_xj array with kernel evals
      // k(xl-xj)
      for(n=0, ib=0; n < blklen; n++){
        for(m=0; m < blklen; m++,ib++){
          blk_xj[ib] = 1.0;


          for(l = 0; l < num_reg_continuous; l++){
            blk_xj[ib] *= xck((matrix_X_continuous[l][blk+n]-matrix_X_continuous[l][blj+m])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
          }

          for(l = 0; l < num_reg_unordered; l++){
            blk_xj[ib] *= xuk(matrix_X_unordered[l][blk+n]==matrix_X_unordered[l][blj+m],
                              lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
          }

          for(l = 0; l < num_reg_ordered; l++){
            blk_xj[ib] *= xok(matrix_X_ordered[l][blk+n],matrix_X_ordered[l][blj+m],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
          }

          // k(yj-yi)
          ts = 1.0;

          for(l = 0; l < num_var_continuous; l++){
            ts *= yck((matrix_Y_continuous[l][blk+n]-matrix_Y_continuous[l][blj+m])/matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
          }


          for(l = 0; l < num_var_unordered; l++){
            ts *= yuk(matrix_Y_unordered[l][blk+n]==matrix_Y_unordered[l][blj+m],lambda[l],num_categories[l]);
          }

          for(l = 0; l < num_var_ordered; l++){
            ts *= yok(matrix_Y_ordered[l][blk+n],matrix_Y_ordered[l][blj+m],lambda[l+num_var_unordered]);
          }

          // accumulate marginals and kernel sums along the way
          if((blk+n >= k) && (blj+m >= j) && (blk+n != blj+m)){
            sum_ker[blk+n] += blk_xj[ib]*ts;
            sum_ker_marginal[blk+n] += blk_xj[ib];
          }

        }
      }

      
      for(i=0; i < num_obs; i+=blklen){
        bli = MIN(num_obs-blklen, i);
        
        // second: fill in k(xl-xi) and k(yj-yi) arrays
        for(n=0, ib=0; n < blklen; n++){
          for(m=0; m < blklen; m++,ib++){
            // k(xl-xi) 
            blk_xi[ib] = 1.0;

            for(l = 0; l < num_reg_continuous; l++){
              blk_xi[ib] *= xck((matrix_X_continuous[l][blk+n]-matrix_X_continuous[l][bli+m])/matrix_bandwidth_reg[l][0])/matrix_bandwidth_reg[l][0];
            }

            for(l = 0; l < num_reg_unordered; l++){
              blk_xi[ib] *= xuk(matrix_X_unordered[l][blk+n]==matrix_X_unordered[l][bli+m],
                                lambda[l+num_var_unordered+num_var_ordered],num_categories[l+num_var_unordered+num_var_ordered]);
            }

            for(l = 0; l < num_reg_ordered; l++){
              blk_xi[ib] *= xok(matrix_X_ordered[l][blk+n],matrix_X_ordered[l][bli+m],lambda[l+num_var_unordered+num_var_ordered+num_reg_unordered]);
            }

            // k(2)(yj-yi)
            blk_yij[ib] = 1.0;

            for(l = 0; l < num_var_continuous; l++){
              blk_yij[ib] *= kernel_convol(KERNEL_den,BANDWIDTH_den,
                                           (matrix_Y_continuous[l][blj+n]-matrix_Y_continuous[l][bli+m])/matrix_bandwidth_var[l][0],matrix_bandwidth_var[l][0],
                                           matrix_bandwidth_var[l][0])/matrix_bandwidth_var[l][0];
            }

            for(l = 0; l < num_var_unordered; l++){
              blk_yij[ib] *= kernel_unordered_convolution(KERNEL_unordered_den, matrix_Y_unordered[l][blj+n],
                                                          matrix_Y_unordered[l][bli+m],lambda[l], num_categories[l], matrix_categorical_vals[l]);
            }

            for(l = 0; l < num_var_ordered; l++){
              blk_yij[ib] *= kernel_ordered_convolution(KERNEL_ordered_den, matrix_Y_ordered[l][blj+n],matrix_Y_ordered[l][bli+m], lambda[l+num_var_unordered], 
                                                        num_categories[l+num_var_unordered], matrix_categorical_vals[l+num_var_unordered]);
            }
          }
        }

        // adjustments for re-centering on num_obs % 64 != 0 data
        ms = i-bli;
        ns = j-blj;
        os = k-blk;

        // third: do partial convolution
        // here there be dragons

        for(o=os, ob=os*blklen; o < blklen; o++, ob+=blklen){ //controls l-indexing
          for(n=ns, ib=ns*blklen; n < blklen; n++, ib+=blklen){ //controls j-indexing
            if(blj+n == blk+o)
              continue;

            ts = 0.0;
            for(m=ms; m < blklen; m++){ //controls i-indexing
              if(bli+m == blk+o)
                continue;
              ts += blk_xi[ob+m]*blk_yij[ib+m];
            }
            sum_ker_convol[blk+o] += blk_xj[ob+n]*ts;
          }
        }
      }
    }
  }


  for(j = 0; j < num_obs; j++){
    if(sum_ker_marginal[j] <= 0.0){
      *cv = FLT_MAX;
      break;
    }
    *cv += (sum_ker_convol[j]/sum_ker_marginal[j]-2.0*sum_ker[j])/sum_ker_marginal[j];
  }

  if (*cv != FLT_MAX)
    *cv /= (double) num_obs;
    
  free(lambda);
  free_mat(matrix_bandwidth_var,num_var_continuous);
  free_mat(matrix_bandwidth_reg,num_reg_continuous);

  free(blk_xi);
  free(blk_xj);
  free(blk_yij);

  free(sum_ker);
  free(sum_ker_convol);
  free(sum_ker_marginal);

  return(0);
}


/*

int np_cuokernelv_loo_mlcv(int KERNEL, int uKERNEL, int oKERNEL,
                           int BANDWIDTH_den,
                           const int num,
                           const int unum, const int onum, const int cnum, 
                           double * const * const xtu,
                           double * const * const xto,
                           double * const * const xtc,
                           double * vector_scale_factor,
                           int * const ncat,                            
                           double * const cv){

  int i,j,k; 

  double (* const ck)(double) = allck[KERNEL];
  double (* const uk)(int, double, int) = alluk[uKERNEL];
  double (* const ok)(double, double, double) = allok[oKERNEL];

  double tmpk, tmps;

  double *lambda;
  double **matrix_bandwidth;

  const double log_FLT_MIN = log(FLT_MIN);

  assert(BANDWIDTH_den == BW_FIXED);

  
  assert((unum <= onum) && (onum <= cnum));

  lambda = alloc_vecd(unum+onum);
  matrix_bandwidth = alloc_matd(num,cnum);

  if(kernel_bandwidth_mean(KERNEL,
                           BANDWIDTH_den,
                           num,
                           num,
                           0,
                           0,
                           0,
                           cnum,
                           unum,
                           onum,
                           vector_scale_factor,
                           xtc,
                           xtc,
                           xtc,
                           xtc,
                           matrix_bandwidth,                         
                           matrix_bandwidth,
                           lambda)==1){
    free(lambda);
    free_mat(matrix_bandwidth,cnum);
    return(1);
  }
  
  *cv = 0.0;

  for(k = 0; k < num; k++){
    tmps = 0.0;

    for(j = 0; j < num; j++){
      if(j == k) continue;

      tmpk = 1.0/num;
      i = 0;
      do{
        switch((i >= unum)+(i >= onum)){
        case 0:tmpk *= uk((xtu[i][j]==xtu[i][k]), lambda[i], ncat[i]);
        case 1:tmpk *= ok(xto[i][j], xto[i][k], lambda[unum+i]);
        case 2:tmpk *= ck((xtc[i][j]-xtc[i][k])/matrix_bandwidth[i][0])/matrix_bandwidth[i][0];
        }
      } while(++i < cnum);
      tmps += tmpk;
    }
    *cv -= (tmps > FLT_MIN) ? log(tmps) : log_FLT_MIN;
  }

  free(lambda);
  free_mat(matrix_bandwidth, cnum);
  return(0);
}


*/
