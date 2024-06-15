#include "gqdouble.h"
#include "operators.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "qlmps/qlmps.h"
#include "singlesiteupdate2.h"
#include "twositeupdate2.h"
#include "myutil.h"
#include "two_site_update_noised_finite_vmps_mpi_impl2.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

#include "params_case.h"

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env(mpi::threading::multiple);
  if (env.thread_level() < mpi::threading::multiple) {
    std::cout << "thread level of env is not right." << std::endl;
    env.abort(-1);
  }
  mpi::communicator world;
  CaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
  size_t N = (1 + Np) * (Lx * Ly);
  float t = params.t, g = params.g, U = params.U, omega = params.omega;

  if (world.rank() == 0) {
    std::cout << "Holstein Model DMRG." << std::endl;
    cout << "System size = (" << Lx << "," << Ly << ")" << endl;
    cout << "The number of electron sites =" << Lx * Ly << endl;
    cout << "The number of phonon pseudosite (per bond) =" << Np << endl;
    cout << "The number of phonon pseudosite (total) =" << Np * Lx * Ly << endl;
    cout << "The total number of sites = " << N << endl;
    cout << "Model parameter: t =" << t << ", g =" << g << ", U =" << U << ",omega=" << omega << endl;
  }
  std::vector<size_t> input_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      input_D_set);

  clock_t startTime, endTime;
  startTime = clock();
  OperatorInitial();
  vector<IndexT> pb_out_set(N);
  vector<long> Tx(N, -1), Ty(N, -1), ElectronSite(Lx * Ly);
  auto iter = ElectronSite.begin();
  // translation along x(for electron) and translation along y(for electron);
  for (int i = 0; i < N; ++i) {
    if (i % (Np + 1) == 0) {
      pb_out_set[i] = pb_outF;
      *iter = i;
      iter++;
    } else pb_out_set[i] = pb_outB;
  }

  SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(pb_out_set);
  MPO<Tensor> mpo(N);
  for (size_t i = 0; i < mpo.size(); i++) {
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
    mpo.LoadTen(i, filename);
  }

  if (world.rank() == 0) {
    cout << "MPO loaded." << endl;
  }
  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
  FiniteMPST mps(sites);

  std::vector<long unsigned int> stat_labs(N, 1);
  int sitenumber_perhole;
  if (params.Numhole > 0) {
    sitenumber_perhole = ElectronSite.size() / params.Numhole;
  } else {
    sitenumber_perhole = 4 * ElectronSite.size();
  }

  int qn_label = 1;
  for (int i = 0; i < ElectronSite.size(); i++) {
    if (i % sitenumber_perhole == sitenumber_perhole / 2) {
      stat_labs[ElectronSite[i]] = 3;
    } else {
      stat_labs[ElectronSite[i]] = qn_label;
      qn_label = 3 - qn_label;
    }
  }

  if (world.rank() == 0) {
    if (params.TotalThreads > 2) {

      qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads - 2);
    } else {

      qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);
    }
  } else {

    qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);
  }

  qlmps::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );

  if (world.rank() == 0) {
    if (IsPathExist(kMpsPath)) {//mps only can be load from file
      if (N == GetNumofMps()) {
        cout << "The number of mps files is consistent with mps size." << endl;
        cout << "Directly use mps from files." << endl;
      } else {
        qlmps::DirectStateInitMps(mps, stat_labs);
        cout << "Initial mps as direct product state." << endl;
        mps.Dump(sweep_params.mps_path, true);
      }
    } else {
      qlmps::DirectStateInitMps(mps, stat_labs);
      cout << "Initial mps as direct product state." << endl;
      mps.Dump(sweep_params.mps_path, true);
    }
  }

  double e0;
  if (!has_bond_dimension_parameter) {
    e0 = qlmps::TwoSiteFiniteVMPS2(mps, mpo, sweep_params, world);
  } else {
    size_t DMRG_time = input_D_set.size();
    std::vector<size_t> MaxLanczIterSet(DMRG_time);
    MaxLanczIterSet.back() = params.MaxLanczIter;
    if (DMRG_time > 1) {
      size_t MaxLanczIterSetSpace;
      MaxLanczIterSet[0] = 3;
      MaxLanczIterSetSpace = (params.MaxLanczIter - 3) / (DMRG_time - 1);
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << ", ";
      for (size_t i = 1; i < DMRG_time - 1; i++) {
        MaxLanczIterSet[i] = MaxLanczIterSet[i - 1] + MaxLanczIterSetSpace;
        std::cout << MaxLanczIterSet[i] << ", ";
      }
      std::cout << MaxLanczIterSet.back() << "]" << std::endl;
    } else {
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << "]" << std::endl;
    }

    for (size_t i = 0; i < DMRG_time; i++) {
      size_t D = input_D_set[i];
      if (world.rank() == 1) {
        std::cout << "D_max = " << D << std::endl;
      }
      sweep_params.Dmax = D;
      sweep_params.Dmin = D;
      e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
    }
  }

  if (world.rank() == 0) {
    std::cout << "E0/site: " << e0 / N << std::endl;
    endTime = clock();
    cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  }
  return 0;
}
