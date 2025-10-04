#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <stdlib.h>
#include <mpi.h>

#include "myutil.h"
#include "my_measure_appendix.h"

#include "qlten/utility/timer.h"

using std::cout;
using std::endl;
using std::vector;
using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
using qlmps::SiteVec;
using qlmps::MeasureOneSiteOp;
using qlten::Timer;
using namespace qlten;
using namespace qlmps;

int main(int argc, char *argv[]) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);
  clock_t startTime, endTime;
  startTime = clock();

  CaseParams params(argv[1]);

  size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
  size_t N = Lx * Ly + (2 * Lx * Ly - Ly) * Np;
  if (GetNumofMps() != N) {
    std::cout << "The number of mps files are inconsistent with mps size!" << std::endl;
    exit(1);
  }

  if (Np != 3) {
    std::cout << "only support Np = 3" << std::endl;
    exit(2);
  }
  OperatorInitial();

  std::vector<IndexT> pb_out_set(N);
  std::vector<size_t> Fsite_set(Lx * Ly);
  std::vector<size_t> Bsite_set(N - Lx * Ly);
  auto iterF = Fsite_set.begin();
  auto iterB = Bsite_set.begin();

  for (size_t i = 0; i < N; ++i) {
    if (IsElectron(i, Ly, Np)) {
      pb_out_set[i] = pb_outF;
      *iterF = i;
      iterF++;
    } else {
      pb_out_set[i] = pb_outB;
      *iterB = i;
      iterB++;
    }
  }
  cout << "The Fermion sites: " << endl;
  Show(Fsite_set);
  cout << '\n';
  cout << "The Boson sites:" << endl;
  Show(Bsite_set);

  SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(pb_out_set);
  FiniteMPST mps(sites);
  qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);

  Timer onesite_timer("measure phonon displacement");

  if (params.Np == 3) {
    MeasuRes<TenElemT> x1, x2, x3, x4, x_res;
    std::vector<double> x1_vec, x2_vec, x3_vec, x4_vec, x_vec;//res

    std::vector<Tensor> op_vec1 = {a, a, 2.0 * adag};
    std::vector<Tensor> op_vec1dag = {adag, adag, 2.0 * a};
    std::vector<Tensor> op_vec2 = {a, adag, sqrt(2) * aadag + sqrt(6) * n_a};
    std::vector<Tensor> op_vec2dag = {adag, a, sqrt(2) * aadag + sqrt(6) * n_a};
    std::vector<Tensor> op_vec3 = {x, n_a, sqrt(3) * aadag + sqrt(7) * n_a};
    std::vector<Tensor> op_vec4 = {x, aadag, aadag + sqrt(5) * n_a};
    switch (rank) {
      case 0:x1 = MeasureOnePhoneOp(mps, op_vec1, Bsite_set, "phonon_x1");
        x1_vec.reserve(x1.size());
        for (size_t i = 0; i < x1.size(); i++) {
          x1_vec.push_back(x1[i].avg);
        }
        x2_vec.resize(x1_vec.size());
        x3_vec.resize(x1_vec.size());
        x4_vec.resize(x1_vec.size());
        MPI_Recv(x2_vec.data(), static_cast<int>(x2_vec.size() * sizeof(double)), MPI_BYTE, 1, 1, comm, MPI_STATUS_IGNORE);
        MPI_Recv(x3_vec.data(), static_cast<int>(x3_vec.size() * sizeof(double)), MPI_BYTE, 2, 2, comm, MPI_STATUS_IGNORE);
        MPI_Recv(x4_vec.data(), static_cast<int>(x4_vec.size() * sizeof(double)), MPI_BYTE, 3, 3, comm, MPI_STATUS_IGNORE);
        break;
      case 1:x2 = MeasureOnePhoneOp(mps, op_vec2, Bsite_set, "phonon_x2");
        x2_vec.reserve(x2.size());
        for (size_t i = 0; i < x2.size(); i++) {
          x2_vec.push_back(x2[i].avg);
        }
        MPI_Send(x2_vec.data(), static_cast<int>(x2_vec.size() * sizeof(double)), MPI_BYTE, 0, 1, comm);
        break;
      case 2:x3 = MeasureOnePhoneOp(mps, op_vec3, Bsite_set, "phonon_x3");
        x3_vec.reserve(x3.size());
        for (size_t i = 0; i < x3.size(); i++) {
          x3_vec.push_back(x3[i].avg);
        }
        MPI_Send(x3_vec.data(), static_cast<int>(x3_vec.size() * sizeof(double)), MPI_BYTE, 0, 2, comm);
        break;
      case 3:x4 = MeasureOnePhoneOp(mps, op_vec4, Bsite_set, "phonon_x4");
        x4_vec.reserve(x4.size());
        for (size_t i = 0; i < x4.size(); i++) {
          x4_vec.push_back(x4[i].avg);
        }
        MPI_Send(x4_vec.data(), static_cast<int>(x4_vec.size() * sizeof(double)), MPI_BYTE, 0, 3, comm);
        break;
      default:break;
    }
    //for matlab plot:
    //x = 2*x1 + 2*x2 + x3 + x4;
    if (rank == 0) {
      x_vec.resize(x1_vec.size());
      for (size_t i = 0; i < x1.size(); i++) {
        x_vec[i] = 2 * x1_vec[i] + 2 * x2_vec[i] + x3_vec[i] + x4_vec[i];
        x_res.push_back(MeasuResElem({i}, x_vec[i]));
      }
    }
    DumpMeasuRes(x_res, "phonon_x");
    std::ofstream ofs("phonon_x", std::ofstream::binary);
    ofs.write((const char *) x_vec.data(), x_vec.size() * sizeof(double));
    ofs << std::endl;
    ofs.close();
  } else {
    std::cout << "not support" << std::endl;
  }
  onesite_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  MPI_Finalize();
  return 0;
}

