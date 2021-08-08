#include<iostream>
#include<cstdlib>
#include "mkl.h"


using std::cout;
using std::endl;
using std::min;

int main()
{ 
  mkl_set_num_threads(10);
  double *A, *B, *C;
  double *A2;
  int m,n,k,i,j;
  double alpha, beta;
  std::cout << "\n This example computes real matrix C=alpha*A*B+beta*C using \n"
            << " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
            << " alpha and beta are double precision scalars\n\n";
  m = 10, k = 10, n = 10;
  cout << " Initializing data for matrix multiplication C=A*B for matrix \n"
       <<  " A("<< m <<"x"<< k<<" and matrix B("<<k<<"x" <<n <<"\n\n";
  alpha = 1.0; beta=1.0;

  cout << " Allocating memory for matrices aligned on 64-byte boundary for better\n"
        << "performace \n\n";
  A = (double *)mkl_malloc(m*k*sizeof(double), 64);
  A2 = (double *)mkl_malloc(m*k*sizeof(double), 64);
  B = (double *)mkl_malloc(k*n*sizeof(double), 64);  
  C = (double *)mkl_malloc(m*n*sizeof(double), 64);
  if( A == nullptr || B ==nullptr || C==nullptr){
    cout << "error: no memory" <<endl;
    mkl_free(A);
    mkl_free(B);
    mkl_free(C);
    return 1;
  }
  MKL_INT group_count = 1;
  const MKL_INT group_size = 20; //when tranfer this paramter into cblas, use its pointer

  // const double* A_const = A;
  // const double* B_const = B;
  // const double* C_const = C;
  cout << " Initializing matrix data\n\n";
  for(i = 0;i<(m*k);i++){
    A[i]=1.0;
    A2[i] = 2.0;
  }
  for(i = 0;i<(k*n);i++){
    B[i]=1;
  }
  for(i = 0;i<(m*n);i++){
    C[i]=0.0;
  }

  decltype(CblasNoTrans)* transa_array = new decltype(CblasNoTrans)[group_size];
  decltype(CblasNoTrans)* transb_array = transa_array;
  MKL_INT* m_array = new MKL_INT[group_size];
  MKL_INT* k_array = new MKL_INT[group_size];
  MKL_INT* n_array = new MKL_INT[group_size];
  const double** a_array = new const double*[group_size];
  const double** b_array = new const double*[group_size];
  double** c_array = new double*[group_size];

  double* alpha_array = new double[group_size];
  double* beta_array = new double[group_size];
  for(size_t i = 0;i<group_size;i++){
    transa_array[i] = CblasNoTrans;
    m_array[i] = m;
    k_array[i] = k;
    n_array[i] = n;
    alpha_array[i] = 1.0;
    beta_array[i] = 1.0;
    a_array[i] = A;
    b_array[i] = B;
    c_array[i] = C;
  }
  // a_array[group_count]



  cout << "computing matrix product add on c batch";
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
  //           m,n,k, alpha, A, k, B,n,beta,C,n);
  cblas_dgemm_batch(CblasRowMajor, transa_array, transb_array,
          m_array,n_array,k_array,alpha_array,
          a_array,k_array,
          b_array,
          n_array,beta_array,c_array,n_array,group_count,
          &group_size);

  printf (" Top left corner of matrix A: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(k,6); j++) {
        printf ("%12.0f", A[j+i*k]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix B: \n");
    for (i=0; i<min(k,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.0f", B[j+i*n]);
      }
      printf ("\n");
    }
    
    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C[j+i*n]);
      }
      printf ("\n");
    }

    printf ("\n Deallocating memory \n\n");
    mkl_free(A);
    mkl_free(B);
    mkl_free(C);

    printf (" Example completed. \n\n");
    return 0;
  return 0;
}

