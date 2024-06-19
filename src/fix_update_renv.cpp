/**
 * Fix single site update, svdinfo==0 collepse
 *
 */

#include "gqdouble.h"
// #include "qlmps/algorithm/lanczos_solver.h"
// #include "qlmps/algorithm/lanczos_solver_impl.h"
#include "qlmps/qlmps.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "myutil.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using namespace qlmps;
using namespace qlten;
using namespace std;

int Parser(const int argc, char *argv[],
           size_t &from,
           size_t &to,
           size_t &thread);

/**
 *
 * @example
 * ./fix_update_renv --thread=24  --from=10 --to=35
 * which uses 24 threads, update from the renv10.qlten (note you require to prepare renv9.qlten and other mps mpo tensors),
 * and update to renv34.qlten (note renv35.qlten will not get);
 *
 * the `from` renv will be generated, and `to` won't.
 *
 */
int main(int argc, char *argv[]) {
  std::cout << "This program used to update right environment tensors by setting the starting/ending sites"
            << std::endl;
  std::cout << "The temp file should named as .temp" << std::endl;
  std::cout << "Note you should give all of the requiring tensors" << std::endl;
  size_t from(0), to(0), thread(0);
  Parser(argc, argv, from, to, thread);

  std::cout << "Argument read: " << std::endl;
  std::cout << "from = " << from << std::endl;
  std::cout << "to = " << to << std::endl;
  std::cout << "thread = " << thread << std::endl;


  qlten::hp_numeric::SetTensorManipulationThreads(thread);

  const size_t N = GetNumofMps();
  using TenT = Tensor;
  const string temp_path = kRuntimeTempPath;
  const size_t Lx(48), Ly(4), Np(3);
  vector<IndexT> pb_out_set(N);
  // translation along x(for electron) and translation along y(for electron);
  for (int i = 0; i < N; ++i) {
    int residue = i % ((2 * Np + 1) * Ly);
    if (residue < (Np + 1) * Ly && residue % (Np + 1) == 0) {
      pb_out_set[i] = pb_outF;
    } else pb_out_set[i] = pb_outB;
  }

  SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(pb_out_set);

  FiniteMPS<TenElemT, U1U1QN> mps(sites);
  MPO<Tensor> mpo(N);
  const std::string kMpoPath = "mpo";
  const std::string kMpoTenBaseName = "mpo_ten";
  for (size_t i = N - from; i > N - to; i--) {
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
    mpo.LoadTen(i, filename);
  }

  cout << "MPO loaded." << endl;
  std::string mps_path = kMpsPath;

  if (from == 0) {
    mps.LoadTen(N - 1, GenMPSTenName(mps_path, N - 1));
    auto mps_trivial_index = mps.back().GetIndexes()[2];
    auto mpo_trivial_index_inv = InverseIndex(mpo.back().GetIndexes()[3]);
    auto mps_trivial_index_inv = InverseIndex(mps_trivial_index);
    TenT renv = TenT({mps_trivial_index, mpo_trivial_index_inv, mps_trivial_index_inv});
    renv({0, 0, 0}) = 1;
    mps.dealloc(N - 1);
    string file = GenEnvTenName("r", 0, temp_path);
    WriteQLTensorTOFile(renv, file);
    from = 1;
  }

  string file = GenEnvTenName("r", from - 1, temp_path);
  TenT renv;
  ifstream renv_file(file);
  renv_file >> renv;

  mps.LoadTen(N - from, GenMPSTenName(mps_path, N - from));
  bool new_code;
  if (mps[N - from].GetIndexes()[2] == renv.GetIndexes()[0]) {
    new_code = true;
  } else if (mps[N - from].GetIndexes()[2] == renv.GetIndexes()[2]) {
    new_code = false;
  } else {
    std::cout << "unexpected, code has bug?" << std::endl;
    exit(1);
  }

  for (size_t i = from; i < to; i++) {
    std::string mps_file_name = GenMPSTenName(mps_path, N - i);
    mps.LoadTen(N - i, mps_file_name);
    auto file = GenEnvTenName("r", i, temp_path);

    if (!new_code) {
      TenT temp1;
      Contract(&mps[N - i], &renv, {{2}, {0}}, &temp1);
      renv = TenT();
      TenT temp2;
      Contract(&temp1, &mpo[N - i], {{1, 2}, {1, 3}}, &temp2);
      auto mps_ten_dag = Dag(mps[N - i]);
      Contract(&temp2, &mps_ten_dag, {{3, 1}, {1, 2}}, &renv);
    } else {
      renv = std::move(UpdateSiteRenvs(renv, mps[N - i], mpo[N - i]));
    }

    WriteQLTensorTOFile(renv, file);
    mps.dealloc(N - i);
  }

  return 0;

}

int Parser(const int argc, char *argv[],
           size_t &from,
           size_t &to,
           size_t &thread) {
  int nOptionIndex = 1;

  string arguement1 = "--from=";
  string arguement2 = "--to=";
  string arguement3 = "--thread=";
  bool from_argument_has(false), to_argument_has(false), thread_argument_has(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement1.size()];
      from = atoi(para_string.c_str());
      from_argument_has = true;
    } else if (strncmp(argv[nOptionIndex], arguement2.c_str(), arguement2.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement2.size()];
      to = atoi(para_string.c_str());
      to_argument_has = true;
    } else if (strncmp(argv[nOptionIndex], arguement3.c_str(), arguement3.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement3.size()];
      thread = atoi(para_string.c_str());
      thread_argument_has = true;
    } else {
      cout << "Options '" << argv[nOptionIndex] << "' not valid. Run '" << argv[0] << "' for details." << endl;
      //   return -1;
    }
    nOptionIndex++;
  }

  if (!from_argument_has) {
    from = 0;
    std::cout << "Note: no from argument, set it as 0 by default." << std::endl;
  }

  if (!to_argument_has) {
    to = 10;
    std::cout << "Note: no to argument, set it as 10 by default." << std::endl;
  }

  if (!thread_argument_has) {
    thread = 24;
    std::cout << "Note: no thread argument, set it as 24 by default." << std::endl;
  }

  return 0;

}