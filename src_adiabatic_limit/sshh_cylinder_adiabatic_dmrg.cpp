/**
 * SSHH model at adiabatic limit
 *
 * data of X_ij are stored in the directory phonon_displacement
 */
#include "gqmps2/gqmps2.h"
#include "./gqdouble.h"
#include "./operators.h"
#include "../src/myutil.h"
#include "./params_case.h"

using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1U1QN>;

int main(int argc, char *argv[]) {
  using namespace std;
  using namespace gqmps2;
  using namespace gqten;
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;
  CaseParams params(argv[1]);

  if (world.rank() == 0 && world.size() > 1 && params.TotalThreads > 2) {
    gqten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads - 2);
    gqten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads - 2);
  } else {
    gqten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads);
    gqten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);
  }

  /******** Model parameter ********/
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = Lx * Ly;
  double t = params.t, alpha = params.alpha, U = params.U, K = params.K;
  double W = 8.0 * t; // band width
  double lambda = alpha * alpha / K / W;
  cout << "System size = (" << Lx << "," << Ly << ")" << endl;
  cout << "The number of electron sites =" << Lx * Ly << endl;
  cout << "Model parameter: t :" << t << ", U :" << U
       << ", alpha :" << alpha << ", K :" << K
       << ", lambda : " << lambda
       << endl;

  /****** DMRG parameter *******/
  gqmps2::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      gqmps2::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );

  /********* phonon displacement data **********/
  double *horizontal_x = new double[N - Ly];  // labels are defined by attached to the left sites
  double *vertical_x = new double[N];   // labels are defined by attached to the top sites
  const std::string phonon_disp_data_dir = "phonon_displacement";
  const std::string horizontal_file = "horizontal_x";
  const std::string vertical_file = "vertical_x";
  if (world.rank() == 0) {
    if (IsPathExist(phonon_disp_data_dir)) {
      //has phonon data, load the data
      //two files, one horizontal_x, has N data,
      //one vertical_x, has N - Ly data
      std::ifstream ifs(phonon_disp_data_dir + horizontal_file, std::ofstream::binary);
      ifs.read((char *) horizontal_x, (N - Ly) * sizeof(double));
      ifs.close();
      ifs.open(phonon_disp_data_dir + vertical_file, std::ofstream::binary);
      ifs.read((char *) vertical_x, N * sizeof(double));
      ifs.close();
    } else {
      // random generate the phonon displacement;
      CreatPath(phonon_disp_data_dir);

      std::default_random_engine random_engine;
      random_engine.seed(std::random_device{}());
      std::uniform_real_distribution<double> u(0, 1);
      for (size_t i = 0; i < N - Ly; i++) {
        horizontal_x[i] = u(random_engine);
      }
      for (size_t i = 0; i < N; i++) {
        vertical_x[i] = u(random_engine);
      }
    }
    // the phonon displacement data are only stored in main processor.
  }

  clock_t startTime, endTime;
  startTime = clock();
  OperatorInitial();
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_out);

  /****** Initialize MPS ******/
  FiniteMPST mps(sites);
  if (world.rank() == 0) {
    if (!IsPathExist(kMpsPath) || !(N == GetNumofMps())) {
      cout << "Initial mps as direct product state." << endl;
      std::vector<size_t> stat_labs(N, 0);
      size_t sitenumber_perhole;
      if (params.Numhole > 0) {
        sitenumber_perhole = N / params.Numhole;
      } else {
        sitenumber_perhole = 4 * N;
      }

      size_t half_filled_qn_label = 1;
      for (size_t i = 0; i < N; i++) {
        if (i % sitenumber_perhole == sitenumber_perhole / 2) {
          stat_labs[i] = 3;   //punch a hole
        } else {
          stat_labs[i] = half_filled_qn_label;
          half_filled_qn_label = 3 - half_filled_qn_label;
        }
      }
      gqmps2::DirectStateInitMps(mps, stat_labs);
      mps.Dump(sweep_params.mps_path, true);
    }
  }

  /******* Set the Measurement Sites ******/
  std::vector<Tensor> hopping1 = {bupcF, bupa},
      hopping2 = {bdnc, Fbdna},
      hopping3 = {-bupaF, bupc},
      hopping4 = {-bdna, Fbdnc};
  std::vector<std::vector<size_t>> horizontal_bond_sites(N - Ly);
  std::vector<std::vector<Tensor>> horizontal_bond_inserts(N - Ly);
  for (size_t i = 0; i < N - Ly; i++) {
    size_t site1 = i, site2 = i + Ly;
    horizontal_bond_sites[i] = {site1, site2};
    horizontal_bond_inserts[i] = std::vector<Tensor>(Ly - 1, f);
  }
  std::vector<std::vector<size_t>> vertical_bond_sites(N);
  std::vector<std::vector<Tensor>> vertical_bond_inserts(N);
  for (size_t i = 0; i < N; i++) {
    size_t y = i % Ly, x = i / Ly;
    if (y < Ly - 1) {
      size_t site1 = i, site2 = i + 1;
      vertical_bond_sites[i] = {site1, site2};
      vertical_bond_inserts[i] = std::vector<Tensor>(0);
    } else { // assume Ly > 2
      size_t site1 = i - Ly + 1, site2 = i;
      vertical_bond_sites[i] = {site1, site2};
      vertical_bond_inserts[i] = std::vector<Tensor>(Ly - 2, f);
    }
  }


  /*******  Calculation Subroutine *******/
  for (size_t iter = 0; iter < params.HamiltonianIters; iter++) {
    // creation the mpo/mro by the newest phonon displacement data
    gqmps2::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);
    for (size_t i = 0; i < N; ++i) {
      mpo_gen.AddTerm(U, Uterm, i);
    }

    //horizontal hopping
    for (size_t i = 0; i < N - Ly; ++i) {
      size_t site1 = i, site2 = i + Ly;
      size_t x = i / Ly;
      mpo_gen.AddTerm(-t + alpha * horizontal_x[i], bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-t + alpha * horizontal_x[i], bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(t - alpha * horizontal_x[i], bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(t - alpha * horizontal_x[i], bdna, site1, Fbdnc, site2, f);
    }
    //vertical hopping
    for (size_t i = 0; i < N; ++i) {
      size_t y = i % Ly, x = i / Ly;
      if (y < Ly - 1) {
        size_t site1 = i, site2 = i + 1;
        mpo_gen.AddTerm(-t + alpha * vertical_x[i], bupcF, site1, bupa, site2);
        mpo_gen.AddTerm(-t + alpha * vertical_x[i], bdnc, site1, Fbdna, site2);
        mpo_gen.AddTerm(t - alpha * vertical_x[i], bupaF, site1, bupc, site2);
        mpo_gen.AddTerm(t - alpha * vertical_x[i], bdna, site1, Fbdnc, site2);
//      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
      } else if (Ly > 2) {
        size_t site1 = i - Ly + 1, site2 = i;
        mpo_gen.AddTerm(-t + alpha * vertical_x[i], bupcF, site1, bupa, site2, f);
        mpo_gen.AddTerm(-t + alpha * vertical_x[i], bdnc, site1, Fbdna, site2, f);
        mpo_gen.AddTerm(t - alpha * vertical_x[i], bupaF, site1, bupc, site2, f);
        mpo_gen.AddTerm(t - alpha * vertical_x[i], bdna, site1, Fbdnc, site2, f);
//      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
      }
    }

    auto mro = mpo_gen.GenMatReprMPO();
    cout << "MRO generated." << endl;

    // dmrg
    double e0 = gqmps2::FiniteDMRG(mps, mro, sweep_params, world);


    //measure hopping and update the phonon displacement
    //TODO: Optimize
    if (world.rank() == 0) {
      mps.Load(sweep_params.mps_path);
      MeasuRes<TenElemT> hopping_data1 = MeasureTwoSiteOp(mps,
                                                          hopping1,
                                                          horizontal_bond_inserts,
                                                          horizontal_bond_sites,
                                                          "horizontal_cupcreationcupannihilation");
      MeasuRes<TenElemT> hopping_data2 = MeasureTwoSiteOp(mps,
                                                          hopping2,
                                                          horizontal_bond_inserts,
                                                          horizontal_bond_sites,
                                                          "horizontal_cdowncreationcdownannihilation");
      MeasuRes<TenElemT> hopping_data3 = MeasureTwoSiteOp(mps,
                                                          hopping3,
                                                          horizontal_bond_inserts,
                                                          horizontal_bond_sites,
                                                          "horizontal_cupannihilationcupcreation");
      MeasuRes<TenElemT> hopping_data4 = MeasureTwoSiteOp(mps,
                                                          hopping4,
                                                          horizontal_bond_inserts,
                                                          horizontal_bond_sites,
                                                          "horizontal_cdownannihilationcdowncreation");

      for (size_t i = 0; i < N - Ly; i++) {
        horizontal_x[i] =
            -alpha / K * (hopping_data1[i].avg + hopping_data2[i].avg + hopping_data3[i].avg + hopping_data4[i].avg);
      }

      hopping_data1 = MeasureTwoSiteOp(mps,
                                       hopping1,
                                       vertical_bond_inserts,
                                       vertical_bond_sites,
                                       "horizontal_cupcreationcupannihilation");
      hopping_data2 = MeasureTwoSiteOp(mps,
                                       hopping2,
                                       vertical_bond_inserts,
                                       vertical_bond_sites,
                                       "horizontal_cdowncreationcdownannihilation");
      hopping_data3 = MeasureTwoSiteOp(mps,
                                       hopping3,
                                       vertical_bond_inserts,
                                       vertical_bond_sites,
                                       "horizontal_cupannihilationcupcreation");
      hopping_data4 = MeasureTwoSiteOp(mps,
                                       hopping4,
                                       vertical_bond_inserts,
                                       vertical_bond_sites,
                                       "horizontal_cdownannihilationcdowncreation");
      for (size_t i = 0; i < N; i++) {
        vertical_x[i] =
            -alpha / K * (hopping_data1[i].avg + hopping_data2[i].avg + hopping_data3[i].avg + hopping_data4[i].avg);
      }
    }

    if (world.rank() == 0) {
      double phonon_energy = 0.0;
      for (size_t i = 0; i < N; i++) {
        phonon_energy += vertical_x[i] * vertical_x[i];
      }
      for (size_t i = 0; i < N - Ly; i++) {
        phonon_energy += horizontal_x[i] * horizontal_x[i];
      }
      phonon_energy *= (K / 2.0);
      std::cout << "Total Energy: " << e0 + phonon_energy << std::endl;
      // print the phonon displacement
      std::cout << "Horizontal Displacement X : \n";
      std::cout << "[";
      for (size_t i = 0; i < N - Ly; i++) {
        std::cout << " " << horizontal_x[i];
      }
      std::cout << "]\n" << std::endl;

      std::cout << "Vertical Displacement X : \n";
      std::cout << "[";
      for (size_t i = 0; i < N; i++) {
        std::cout << " " << vertical_x[i];
      }
      std::cout << "]\n" << std::endl;

      endTime = clock();
      cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    }
  }


  /****** Dump phonon displacement data *******/
  std::ofstream ofs(phonon_disp_data_dir + horizontal_file, std::ofstream::binary);
  ofs.write((char *) horizontal_x, (N - Ly) * sizeof(double));
  ofs.close();
  ofs.open(phonon_disp_data_dir + vertical_file, std::ofstream::binary);
  ofs.write((char *) vertical_x, N * sizeof(double));
  ofs.close();

  delete[] horizontal_x;
  delete[] vertical_x;
  return 0;
}
