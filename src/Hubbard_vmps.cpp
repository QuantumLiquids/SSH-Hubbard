#include <iostream>
#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "params_case.h"
#include "gqdouble.h"
#include "myutil.h"
#include "boost/mpi.hpp"
using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;

  CaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = Lx * Ly;

  if (world.rank() == 0) {
    cout << "System size = (" << Lx << "," << Ly << ")" << endl;
    cout << "The number of electron sites =" << Lx * Ly << endl;
    cout << "The total number of sites = " << N << endl;
    float t = params.t, g = params.g, U = params.U, omega = params.omega;
    cout << "Model parameter: t =" << t << ", g =" << g << ", U =" << U << ",omega=" << omega << endl;
  }

  clock_t startTime, endTime;
  startTime = clock();

  qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);

  qlmps::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter)
  );

  std::vector<size_t> input_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      input_D_set);

  size_t DMRG_time = input_D_set.size();
  std::vector<size_t> MaxLanczIterSet(DMRG_time);
  if (has_bond_dimension_parameter) {
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
  }

  bool noise_valid(false);
  for (size_t i = 0; i < params.noise.size(); i++) {
    if (params.noise[i] != 0) {
      noise_valid = true;
      break;
    }
  }

  double e0(0.0); //energy

  //const size_t Ly = params.Ly;
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_outF);
  qlmps::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  std::string kMpoPath = "mpo";
  const std::string kMpoTenBaseName = "mpo_ten";

  qlmps::MPO<Tensor> mpo(N);
  if (IsPathExist(kMpoPath)) {
    for (size_t i = 0; i < mpo.size(); i++) {
      std::string filename = kMpoPath + "/" +
          kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
      mpo.LoadTen(i, filename);
    }

    cout << "MPO loaded." << endl;
  } else {
    cout << "No mpo directory. exiting" << std::endl;
    exit(0);
  }

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
  FiniteMPST mps(sites);

  std::vector<size_t> stat_labs(N);
  size_t site_number_per_hole;
  if (params.Numhole > 0) {
    site_number_per_hole = N / params.Numhole;
  } else {
    site_number_per_hole = N + 99;
  }

  size_t sz_label = 0;

  for (size_t i = 0; i < N; ++i) {
    if (i % site_number_per_hole == site_number_per_hole - 1) {
      stat_labs[i] = 3;
    } else {
      stat_labs[i] = sz_label % 2 + 1;
      sz_label++;
    }
  }
  if (world.rank() == 0) {
    if (IsPathExist(kMpsPath)) {
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

  world.barrier();
  if (!has_bond_dimension_parameter) {
    e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
  } else {
    for (size_t i = 0; i < DMRG_time; i++) {
      size_t D = input_D_set[i];
      if (world.rank() == 1) {
        std::cout << "D_max = " << D << std::endl;
      }
      qlmps::FiniteVMPSSweepParams sweep_params(
          params.Sweeps,
          D, D, params.CutOff,
          qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i])
      );
      e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
    }
  }

  if (world.rank() == 0) {
    std::cout << "E0/site: " << e0 / N << std::endl;
  }
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
