/*
    phonon_renormalization.cpp
    EVD phonon density matrix
*/
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <stdlib.h>

#include "myutil.h"
#include "my_measure.h"

#include "qlten/utility/timer.h"

using std::cout;
using std::endl;
using std::vector;
using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
using qlmps::SiteVec;
using qlmps::MeasureOneSiteOp;
using qlten::Timer;
using namespace qlten;

int main(int argc, char *argv[]) {
  clock_t startTime, endTime;
  startTime = clock();

  CaseParams params(argv[1]);

  size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
  size_t N = Lx * Ly + (2 * Lx * Ly - Ly) * Np;
  if (GetNumofMps() != N) {
    std::cout << "The number of mps files are inconsistent with mps size!" << std::endl;
    exit(1);
  }

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
  qlten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads);
  qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);

  size_t left_boundary = qlmps::FindLeftBoundary(mps);
  if (left_boundary + 1 >= Np) {//load all the phonon

  } else {
    for (size_t i = left_boundary + 2; i <= Np; i++) {
      mps.LoadTen(i, qlmps::kMpsPath);
    }
  }

  for (size_t i = left_boundary + 1; i > 1; i--) {
    mps.RightCanonicalizeTen(i);
  }

  Tensor temp = mps[1];
  for (size_t i = 2; i <= Np; i++) {
    Tensor temp1;
    size_t last_idx = temp.Rank() - 1;
    qlten::Contract(&temp, mps(i), {{last_idx}, {0}}, &temp1);
    temp = std::move(temp1);
  }

  Tensor temp_dag = qlten::Dag(temp);
  Tensor density_matrix;
  size_t last_idx = temp.Rank() - 1;
  Contract(&temp, &temp_dag, {{0, last_idx}, {0, last_idx}}, &density_matrix);

  Tensor u, s, vt;
  qlmps::mock_qlten::SVD(&density_matrix, Np, qn0, &u, &s, &vt);

  std::cout << "[";
  for (size_t i = 0; i < s.GetShape()[0]; i++) {
    if (i < s.GetShape()[0] - 1) {
      std::cout << s({i, i}) << ", ";
    } else {
      std::cout << s({i, i}) << "] " << std::endl;
    }
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}

