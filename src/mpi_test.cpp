//
// Created by Hao-Xin on 2021/11/15.
//

#include "mpi.h"
#include <cstdlib>
#include <iostream>
#include <climits>


int main(int argc, char* argv[]) {
  std::cout << "INT_MAX = " << INT_MAX << std::endl;
  int *provided(nullptr);
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, provided);
  const size_t num_data = 1e8 * atoi(argv[1]) ;
  double* B = new double[num_data];
  double start, end;
  std::cout <<"data size = " << num_data << std::endl;
  std::cout <<"start to broadcast" << std::endl;

  start = MPI_Wtime();
  MPI_Bcast((void*)B, num_data, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  end = MPI_Wtime();
  printf("Runtime = %f\n", end-start);
  delete[] B;
  MPI_Finalize();
  return 0;
}