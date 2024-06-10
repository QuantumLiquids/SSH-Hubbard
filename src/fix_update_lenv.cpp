/**
 * Fix left enviroment tensors
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
 * which uses 24 threads, update from the lenv10.qlten (note you require to prepare lenv9.qlten and other mps mpo tensors),
 * and update to lenv34.qlten (and lenv35.qlten will not get);
 *
 * the `from` lenv will be generated, and `to` won't.
 *
 */
int main(int argc, char *argv[]) {
  std::cout << "This program used to update left environment tensors by setting the starting/ending sites" << std::endl;
  std::cout << "The temp file should named as .temp" << std::endl;
  std::cout << "Note you should give all of the requiring tensors" << std::endl;
  size_t from(0), to(0), thread(0);
  Parser(argc, argv, from, to, thread);

  std::cout << "Argument read: " << std::endl;
  std::cout << "from = " << from << std::endl;
  std::cout << "to = " << to << std::endl;
  std::cout << "thread = " << thread << std::endl;

  qlten::hp_numeric::SetTensorTransposeNumThreads(thread);
  qlten::hp_numeric::SetTensorManipulationThreads(thread);

  const size_t N = GetNumofMps();
  using TenT = Tensor;
  const string temp_path = kRuntimeTempPath;
  const size_t Ly(4), Lx(N / Ly), Np(3);
  std::cout << "Note: by default we set Ly = 4, and Np = 3" << "\n";
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
  for (size_t i = from - 1; i < to - 1; i++) {
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
    mpo.LoadTen(i, filename);
  }

  cout << "MPO loaded." << endl;

  std::string mps_path = kMpsPath;
  if (from == 0) {
    mps.LoadTen(0, GenMPSTenName(mps_path, 0));
    auto mps_trivial_index = mps.front().GetIndexes()[0];
    auto mpo_trivial_index_inv = InverseIndex(mpo.front().GetIndexes()[0]);
    auto mps_trivial_index_inv = InverseIndex(mps_trivial_index);
    TenT lenv = TenT({mps_trivial_index, mpo_trivial_index_inv, mps_trivial_index_inv});
    lenv({0, 0, 0}) = 1;
    mps.dealloc(0);
    std::string file = GenEnvTenName("l", 0, temp_path);
    WriteQLTensorTOFile(lenv, file);
    from = 1;
  }

  string file = GenEnvTenName("l", from - 1, temp_path);
  TenT lenv;
  ifstream lenv_file(file);
  lenv_file >> lenv;

  mps.LoadTen(from - 1, GenMPSTenName(mps_path, from - 1));
  bool new_code;
  if (mps[from - 1].GetIndexes()[0] == lenv.GetIndexes()[0]) {
    new_code = true;
  } else if (mps[from - 1].GetIndexes()[0] == lenv.GetIndexes()[2]) {
    new_code = false;
  } else {
    std::cout << "unexpected, code has bug?" << std::endl;
    exit(1);
  }

  for (size_t i = from; i < to; i++) {
    mps.LoadTen(i - 1, GenMPSTenName(mps_path, i - 1));
    auto file = GenEnvTenName("l", i, temp_path);
    if (!new_code) {
      TenT temp1;
      Contract(&mps[i - 1], &lenv, {{0}, {0}}, &temp1);
      lenv = TenT();
      TenT temp2;
      Contract(&temp1, &mpo[i - 1], {{0, 2}, {1, 0}}, &temp2);
      auto mps_ten_dag = Dag(mps[i - 1]);
      Contract(&temp2, &mps_ten_dag, {{1, 2}, {0, 1}}, &lenv);
    } else {
      lenv = std::move(UpdateSiteLenvs(lenv, mps[i - 1], mpo[i - 1]));
    }

    WriteQLTensorTOFile(lenv, file);
    mps.dealloc(i - 1);
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