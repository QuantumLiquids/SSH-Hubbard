#include "gqdouble.h"
#include "operators.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "qlmps/qlmps.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

#include "params_case.h"

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);

  const size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
  const size_t N = Lx * Ly + (2 * Lx * Ly - Ly) * Np;
  cout << "System size = (" << Lx << "," << Ly << ")" << endl;
  cout << "The number of electron sites =" << Lx * Ly << endl;
  cout << "The number of phonon pseudosite (per bond) =" << Np << endl;
  cout << "The number of phonon pseudosite (total) =" << (2 * Lx * Ly - Ly) * Np << endl;
  cout << "The total number of sites = " << N << endl;
  float t = params.t, g = params.g, U = params.U, omega = params.omega;
  cout << "Model parameter: t =" << t << ", g =" << g << ", U =" << U << ",omega=" << omega << endl;

  clock_t startTime, endTime;
  startTime = clock();

  OperatorInitial();

  vector<IndexT> pb_out_set(N);
  vector<long> Tx(N, -1), Ty(N, -1), ElectronSite(Lx * Ly);
  // translation along x(for electron) and translation along y(for electron);
  for (size_t i = 0; i < N; ++i) {
    size_t residue = i % ((2 * Np + 1) * Ly);
    if (residue < (Np + 1) * Ly && residue % (Np + 1) == 0) {
      pb_out_set[i] = pb_outF;
    } else pb_out_set[i] = pb_outB;
  }

  SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(pb_out_set);
  qlmps::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  auto iter = ElectronSite.begin();
  for (unsigned i = 0; i < N; ++i) {
    unsigned residue = i % ((2 * Np + 1) * Ly);
    if (residue < (Np + 1) * Ly && residue % (Np + 1) == 0) {
      mpo_gen.AddTerm(U, Uterm, i);
      std::cout << "add site" << i << "Hubbard U term" << endl;
      *iter = i;
      iter++;
      if (i < N - ((Np + 1) * Ly)) Tx[i] = i + (2 * Np + 1) * Ly;
      if (residue != (Np + 1) * (Ly - 1)) Ty[i] = i + Np + 1;
      else Ty[i] = i - (Np + 1) * (Ly - 1);
    } else if (residue < (Np + 1) * Ly) {
      unsigned residue2 = residue % (Np + 1);
      mpo_gen.AddTerm(omega, (TenElemT) (pow(2, residue2 - 1)) * n_a, i);
      cout << "add site" << i << "phonon potential term (phonon number = "
           << pow(2, residue2 - 1) << ")" << endl;
    } else {
      residue = residue - (Np + 1) * Ly;
      unsigned residue2 = residue % Np;
      mpo_gen.AddTerm(omega, (TenElemT) (pow(2, residue2)) * n_a, i);
      cout << "add site" << i << "phonon potential term (phonon number = "
           << pow(2, residue2) << ")" << endl;
    }
  }

  for (unsigned i = 0; i < N; i++) {
    if (Tx[i] != -1) {
      size_t site1(i), site2(Tx[i]);
      vector<size_t> inst_idxs(Ly - 1);
      auto iter = find(ElectronSite.cbegin(), ElectronSite.cend(), site1);
      for (size_t j = 1; j < Ly; j++) {
        inst_idxs[j - 1] = *(iter + j);
      }
      std::cout << "add site (" << site1 << "," << site2 << ")  hopping  term" << endl;
      mpo_gen.AddTerm(-t, bupcF, site1, bupa, site2, f, inst_idxs);
      mpo_gen.AddTerm(-t, bdnc, site1, Fbdna, site2, f, inst_idxs);
      mpo_gen.AddTerm(t, bupaF, site1, bupc, site2, f, inst_idxs);
      mpo_gen.AddTerm(t, bdna, site1, Fbdnc, site2, f, inst_idxs);
    }
    if (Ty[i] != -1) {
      if (Ty[i] > i) {
        unsigned site1(i), site2(Ty[i]);
        cout << "add site (" << site1 << "," << site2 << ")  hopping  term" << endl;
        mpo_gen.AddTerm(-t, bupcF, site1, bupa, site2);
        mpo_gen.AddTerm(-t, bdnc, site1, Fbdna, site2);
        mpo_gen.AddTerm(t, bupaF, site1, bupc, site2);
        mpo_gen.AddTerm(t, bdna, site1, Fbdnc, site2);
      } else {
        size_t site1(Ty[i]), site2(i);
        vector<size_t> inst_idxs(Ly - 2);
        for (size_t j = 1; j < Ly - 1; j++) {
          inst_idxs[j - 1] = site1 + j * (Np + 1);
        }
        std::cout << "add site (" << site1 << "," << site2 << ")  hopping  term" << endl;
        mpo_gen.AddTerm(-t, bupcF, site1, bupa, site2, f, inst_idxs);
        mpo_gen.AddTerm(-t, bdnc, site1, Fbdna, site2, f, inst_idxs);
        mpo_gen.AddTerm(t, bupaF, site1, bupc, site2, f, inst_idxs);
        mpo_gen.AddTerm(t, bdna, site1, Fbdnc, site2, f, inst_idxs);
      }
    }
  }

  for (size_t i = 0; i < N; i++) {
    //electron phonon interaction, in y-direction bond
    if (Ty[i] != -1 && Ty[i] == (int) (i + Np + 1)) {
      vector<size_t> op_sites(Np + 2);
      for (size_t j = 0; j < Np + 2; j++) {
        op_sites[j] = i + j;
      }
      switch (Np) {
        case 1:mpo_gen.AddTerm(g, {bupcF, x, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, x, Fbdna}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, x, Fbdnc}, op_sites);
          cout << "sites ( " << i << "," << i + 1 << "," << i + 2 << ") electron-phonon interaction term" << endl;
          break;
        case 2: {
          auto op1 = (sqrt(3) - 1) * n_a + idB + sqrt(2) * a;
          auto op2 = sqrt(2) * (adag + (-a));
          mpo_gen.AddTerm(g, {bupcF, x, op1, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, a, op2, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, x, op1, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, a, op2, Fbdna}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, op1, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, a, op2, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, x, op1, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, a, op2, Fbdnc}, op_sites);
        }
          break;
        case 3: {
          auto op1 = sqrt(2) * aadag + sqrt(6) * n_a;
          auto op2 = sqrt(3) * aadag + sqrt(7) * n_a;
          auto op3 = aadag + sqrt(5) * n_a;
          mpo_gen.AddTerm(g, {bupcF, a, a, 2.0 * adag, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, adag, adag, 2.0 * a, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, a, adag, op1, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, adag, a, op1, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, x, n_a, op2, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, x, aadag, op3, bupa}, op_sites);

          mpo_gen.AddTerm(g, {bdnc, a, a, 2.0 * adag, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, adag, adag, 2.0 * a, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, a, adag, op1, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, adag, a, op1, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, x, n_a, op2, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, x, aadag, op3, Fbdna}, op_sites);

          mpo_gen.AddTerm(-g, {bupaF, a, a, 2.0 * adag, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, adag, adag, 2.0 * a, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, a, adag, op1, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, adag, a, op1, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, n_a, op2, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, aadag, op3, bupc}, op_sites);

          mpo_gen.AddTerm(-g, {bdna, a, a, 2.0 * adag, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, adag, adag, 2.0 * a, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, a, adag, op1, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, adag, a, op1, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, x, n_a, op2, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, x, aadag, op3, Fbdnc}, op_sites);
        }
          break;
        case 4: {
          auto m = aadag;
          auto n = n_a;
          auto &P0 = m;
          auto &P1 = n;
          //1
          mpo_gen.AddTerm(g, {bupcF, x, m + sqrt(3) * n, m, m, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, x, sqrt(5) * m + sqrt(7) * n, n, m, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, x, sqrt(9) * m + sqrt(11) * n, m, n, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, x, sqrt(13) * m + sqrt(15) * n, n, n, bupa}, op_sites);

          mpo_gen.AddTerm(g, {bupcF, a, adag, sqrt(2) * m + sqrt(6) * n, m, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, adag, a, sqrt(2) * m + sqrt(6) * n, m, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, a, adag, sqrt(10) * m + sqrt(14) * n, n, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, adag, a, sqrt(10) * m + sqrt(14) * n, n, bupa}, op_sites);

          mpo_gen.AddTerm(g, {bupcF, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, bupa}, op_sites);

          mpo_gen.AddTerm(sqrt(8) * g, {bupcF, a, a, a, adag, bupa}, op_sites);
          mpo_gen.AddTerm(sqrt(8) * g, {bupcF, adag, adag, adag, a, bupa}, op_sites);

          //2
          mpo_gen.AddTerm(g, {bdnc, x, m + sqrt(3) * n, m, m, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, x, sqrt(5) * m + sqrt(7) * n, n, m, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, x, sqrt(9) * m + sqrt(11) * n, m, n, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, x, sqrt(13) * m + sqrt(15) * n, n, n, Fbdna}, op_sites);

          mpo_gen.AddTerm(g, {bdnc, a, adag, sqrt(2) * m + sqrt(6) * n, m, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, adag, a, sqrt(2) * m + sqrt(6) * n, m, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, a, adag, sqrt(10) * m + sqrt(14) * n, n, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, adag, a, sqrt(10) * m + sqrt(14) * n, n, Fbdna}, op_sites);

          mpo_gen.AddTerm(g, {bdnc, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, Fbdna}, op_sites);
          mpo_gen.AddTerm(g, {bdnc, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, Fbdna}, op_sites);

          mpo_gen.AddTerm(sqrt(8) * g, {bdnc, a, a, a, adag, Fbdna}, op_sites);
          mpo_gen.AddTerm(sqrt(8) * g, {bdnc, adag, adag, adag, a, Fbdna}, op_sites);

          //3
          mpo_gen.AddTerm(-g, {bupaF, x, m + sqrt(3) * n, m, m, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, sqrt(5) * m + sqrt(7) * n, n, m, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, sqrt(9) * m + sqrt(11) * n, m, n, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, sqrt(13) * m + sqrt(15) * n, n, n, bupc}, op_sites);

          mpo_gen.AddTerm(-g, {bupaF, a, adag, sqrt(2) * m + sqrt(6) * n, m, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, adag, a, sqrt(2) * m + sqrt(6) * n, m, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, a, adag, sqrt(10) * m + sqrt(14) * n, n, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, adag, a, sqrt(10) * m + sqrt(14) * n, n, bupc}, op_sites);

          mpo_gen.AddTerm(-g, {bupaF, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, bupc}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, bupc}, op_sites);

          mpo_gen.AddTerm(-sqrt(8) * g, {bupaF, a, a, a, adag, bupc}, op_sites);
          mpo_gen.AddTerm(-sqrt(8) * g, {bupaF, adag, adag, adag, a, bupc}, op_sites);

          //4
          mpo_gen.AddTerm(-g, {bdna, x, m + sqrt(3) * n, m, m, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, x, sqrt(5) * m + sqrt(7) * n, n, m, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, x, sqrt(9) * m + sqrt(11) * n, m, n, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, x, sqrt(13) * m + sqrt(15) * n, n, n, Fbdnc}, op_sites);

          mpo_gen.AddTerm(-g, {bdna, a, adag, sqrt(2) * m + sqrt(6) * n, m, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, adag, a, sqrt(2) * m + sqrt(6) * n, m, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, a, adag, sqrt(10) * m + sqrt(14) * n, n, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, adag, a, sqrt(10) * m + sqrt(14) * n, n, Fbdnc}, op_sites);

          mpo_gen.AddTerm(-g, {bdna, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-g, {bdna, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, Fbdnc}, op_sites);

          mpo_gen.AddTerm(-sqrt(8) * g, {bdna, a, a, a, adag, Fbdnc}, op_sites);
          mpo_gen.AddTerm(-sqrt(8) * g, {bdna, adag, adag, adag, a, Fbdnc}, op_sites);
        }
          break;
        default:cout << "This progress does not support for Np > 4 cases" << endl;
          exit(0);
      }
    } else if (Ty[i] != -1) {
      vector<size_t> op_sites(Np + 2);
      op_sites[0] = Ty[i];
      op_sites[1] = i;
      for (size_t j = 0; j < Np; j++) {
        op_sites[j + 2] = i + j + 1;
      }
      vector<Tensor> inst_ops(Np + 1, f);
      vector<vector<size_t>> inst_ops_idxs_set(Np + 1);
      for (size_t j = 1; j < Ly - 1; j++) {
        inst_ops_idxs_set[0].push_back(Ty[i] + j * (Np + 1));
      }
      switch (Np) {
        case 1:mpo_gen.AddTerm(g, {bupcF, bupa, x}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, x}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, x}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, x}, op_sites, inst_ops, inst_ops_idxs_set);
          cout << "sites ( " << Ty[i] << "," << i << "," << i + 1 << ") electron-phonon interaction term" << endl;
          break;
        case 2: {
          auto op1 = (sqrt(3) - 1) * n_a + idB + sqrt(2) * a;
          auto op2 = sqrt(2) * (adag + (-a));
          mpo_gen.AddTerm(g, {bupcF, bupa, x, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, x, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, x, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, x, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
        }
          break;
        case 3: {
          auto op1 = sqrt(2) * aadag + sqrt(6) * n_a;
          auto op2 = sqrt(3) * aadag + sqrt(7) * n_a;
          auto op3 = aadag + sqrt(5) * n_a;
          mpo_gen.AddTerm(g, {bupcF, bupa, a, a, 2.0 * adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, adag, adag, 2.0 * a}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, a, adag, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, adag, a, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, x, n_a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, x, aadag, op3}, op_sites, inst_ops, inst_ops_idxs_set);

          mpo_gen.AddTerm(g, {bdnc, Fbdna, a, a, 2.0 * adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, adag, adag, 2.0 * a}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, a, adag, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, adag, a, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, x, n_a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, x, aadag, op3}, op_sites, inst_ops, inst_ops_idxs_set);

          mpo_gen.AddTerm(-g, {bupaF, bupc, a, a, 2.0 * adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, adag, adag, 2.0 * a}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, a, adag, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, adag, a, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, x, n_a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, x, aadag, op3}, op_sites, inst_ops, inst_ops_idxs_set);

          mpo_gen.AddTerm(-g, {bdna, Fbdnc, a, a, 2.0 * adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, adag, adag, 2.0 * a}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, a, adag, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, adag, a, op1}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, x, n_a, op2}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, x, aadag, op3}, op_sites, inst_ops, inst_ops_idxs_set);
        }
          break;
        case 4: {
          auto m = aadag;
          auto n = n_a;
          auto &P0 = m;
          auto &P1 = n;

          //1
          mpo_gen.AddTerm(g, {bupcF, bupa, x, m + sqrt(3) * n, m, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, x, sqrt(5) * m + sqrt(7) * n, n, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, bupa, x, sqrt(9) * m + sqrt(11) * n, m, n}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, bupa, x, sqrt(13) * m + sqrt(15) * n, n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bupcF, bupa, a, adag, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, bupa, adag, a, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, bupa, a, adag, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, bupa, adag, a, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bupcF, bupa, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, bupa, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(sqrt(8) * g, {bupcF, bupa, a, a, a, adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(sqrt(8) * g, {bupcF, bupa, adag, adag, adag, a}, op_sites, inst_ops, inst_ops_idxs_set);

          //2
          mpo_gen.AddTerm(g, {bdnc, Fbdna, x, m + sqrt(3) * n, m, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, x, sqrt(5) * m + sqrt(7) * n, n, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, Fbdna, x, sqrt(9) * m + sqrt(11) * n, m, n}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, Fbdna, x, sqrt(13) * m + sqrt(15) * n, n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bdnc, Fbdna, a, adag, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, Fbdna, adag, a, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, Fbdna, a, adag, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, Fbdna, adag, a, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bdnc, Fbdna, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, Fbdna, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(sqrt(8) * g, {bdnc, Fbdna, a, a, a, adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(sqrt(8) * g, {bdnc, Fbdna, adag, adag, adag, a}, op_sites, inst_ops, inst_ops_idxs_set);

          //3
          mpo_gen.AddTerm(-g, {bupaF, bupc, x, m + sqrt(3) * n, m, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, bupc, x, sqrt(5) * m + sqrt(7) * n, n, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, x, sqrt(9) * m + sqrt(11) * n, m, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, x, sqrt(13) * m + sqrt(15) * n, n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, a, adag, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, adag, a, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, a, adag, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, adag, a, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, bupc, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-sqrt(8) * g, {bupaF, bupc, a, a, a, adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-sqrt(8) * g, {bupaF, bupc, adag, adag, adag, a}, op_sites, inst_ops, inst_ops_idxs_set);

          //4
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, x, m + sqrt(3) * n, m, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, Fbdnc, x, sqrt(5) * m + sqrt(7) * n, n, m}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, x, sqrt(9) * m + sqrt(11) * n, m, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, x, sqrt(13) * m + sqrt(15) * n, n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, a, adag, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, adag, a, sqrt(2) * m + sqrt(6) * n, m},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, a, adag, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, adag, a, sqrt(10) * m + sqrt(14) * n, n},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, Fbdnc, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-sqrt(8) * g, {bdna, Fbdnc, a, a, a, adag}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-sqrt(8) * g, {bdna, Fbdnc, adag, adag, adag, a}, op_sites, inst_ops, inst_ops_idxs_set);
        }
          break;
        default:cout << "This progress does not support for Np > 4 cases" << endl;
          exit(0);
      }
    }
    // electron phonon interaction, in x-direction bond
    if (Tx[i] != -1) {
      size_t residue = i % ((2 * Np + 1) * Ly);
      size_t y = residue / (Np + 1);
      size_t x_coor = i / ((2 * Np + 1) * Ly);
      size_t phonon1 = x_coor * ((2 * Np + 1) * Ly) + Ly * (Np + 1) + Np * y;
      vector<size_t> op_sites(Np + 2);
      op_sites[0] = i;
      op_sites.back() = Tx[i];
      for (size_t j = 0; j < Np; j++) {
        op_sites[j + 1] = phonon1 + j;
      }
      vector<Tensor> inst_ops(Np + 1, f);
      vector<vector<size_t>> inst_ops_idxs_set(Np + 1);
      auto iter = find(ElectronSite.cbegin(), ElectronSite.cend(), i);
      for (size_t j = 1; j < Ly; j++) {
        size_t inst_idxs = *(iter + j);
        if (inst_idxs < phonon1) {
          inst_ops_idxs_set[0].push_back(inst_idxs);
        } else {
          inst_ops_idxs_set.back().push_back(inst_idxs);
        }

      }
      switch (Np) {
        case 1:mpo_gen.AddTerm(g, {bupcF, x, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, x, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, x, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, x, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          cout << "sites ( " << i << "," << phonon1 << "," << Tx[i] << ") electron-phonon interaction term" << endl;
          break;
        case 2: {
          auto op1 = (sqrt(3) - 1) * n_a + idB + sqrt(2) * a;
          auto op2 = sqrt(2) * (adag + (-a));
          mpo_gen.AddTerm(g, {bupcF, x, op1, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, a, op2, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, x, op1, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, a, op2, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, x, op1, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, a, op2, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, x, op1, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, a, op2, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
        }
          break;
        case 3: {
          auto op1 = sqrt(2) * aadag + sqrt(6) * n_a;
          auto op2 = sqrt(3) * aadag + sqrt(7) * n_a;
          auto op3 = aadag + sqrt(5) * n_a;
          mpo_gen.AddTerm(g, {bupcF, a, a, 2.0 * adag, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, adag, adag, 2.0 * a, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, a, adag, op1, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, adag, a, op1, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, x, n_a, op2, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, x, aadag, op3, bupa}, op_sites, inst_ops, inst_ops_idxs_set);

          mpo_gen.AddTerm(g, {bdnc, a, a, 2.0 * adag, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, adag, adag, 2.0 * a, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, a, adag, op1, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, adag, a, op1, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, x, n_a, op2, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, x, aadag, op3, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);

          mpo_gen.AddTerm(-g, {bupaF, a, a, 2.0 * adag, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, adag, adag, 2.0 * a, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, a, adag, op1, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, adag, a, op1, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, x, n_a, op2, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, x, aadag, op3, bupc}, op_sites, inst_ops, inst_ops_idxs_set);

          mpo_gen.AddTerm(-g, {bdna, a, a, 2.0 * adag, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, adag, adag, 2.0 * a, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, a, adag, op1, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, adag, a, op1, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, x, n_a, op2, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, x, aadag, op3, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
        }
          break;
        case 4: {
          auto m = aadag;
          auto n = n_a;
          auto &P0 = m;
          auto &P1 = n;

          //1
          mpo_gen.AddTerm(g, {bupcF, x, m + sqrt(3) * n, m, m, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, x, sqrt(5) * m + sqrt(7) * n, n, m, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bupcF, x, sqrt(9) * m + sqrt(11) * n, m, n, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, x, sqrt(13) * m + sqrt(15) * n, n, n, bupa},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bupcF, a, adag, sqrt(2) * m + sqrt(6) * n, m, bupa},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, adag, a, sqrt(2) * m + sqrt(6) * n, m, bupa},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, a, adag, sqrt(10) * m + sqrt(14) * n, n, bupa},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, adag, a, sqrt(10) * m + sqrt(14) * n, n, bupa},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bupcF, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, bupa},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bupcF, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, bupa},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(sqrt(8) * g, {bupcF, a, a, a, adag, bupa}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(sqrt(8) * g, {bupcF, adag, adag, adag, a, bupa}, op_sites, inst_ops, inst_ops_idxs_set);

          //2
          mpo_gen.AddTerm(g, {bdnc, x, m + sqrt(3) * n, m, m, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, x, sqrt(5) * m + sqrt(7) * n, n, m, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g, {bdnc, x, sqrt(9) * m + sqrt(11) * n, m, n, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, x, sqrt(13) * m + sqrt(15) * n, n, n, Fbdna},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bdnc, a, adag, sqrt(2) * m + sqrt(6) * n, m, Fbdna},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, adag, a, sqrt(2) * m + sqrt(6) * n, m, Fbdna},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, a, adag, sqrt(10) * m + sqrt(14) * n, n, Fbdna},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, adag, a, sqrt(10) * m + sqrt(14) * n, n, Fbdna},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(g,
                          {bdnc, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, Fbdna},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(g,
                          {bdnc, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, Fbdna},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(sqrt(8) * g, {bdnc, a, a, a, adag, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(sqrt(8) * g, {bdnc, adag, adag, adag, a, Fbdna}, op_sites, inst_ops, inst_ops_idxs_set);

          //3
          mpo_gen.AddTerm(-g, {bupaF, x, m + sqrt(3) * n, m, m, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bupaF, x, sqrt(5) * m + sqrt(7) * n, n, m, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, x, sqrt(9) * m + sqrt(11) * n, m, n, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, x, sqrt(13) * m + sqrt(15) * n, n, n, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bupaF, a, adag, sqrt(2) * m + sqrt(6) * n, m, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, adag, a, sqrt(2) * m + sqrt(6) * n, m, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, a, adag, sqrt(10) * m + sqrt(14) * n, n, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, adag, a, sqrt(10) * m + sqrt(14) * n, n, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bupaF, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bupaF, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, bupc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-sqrt(8) * g, {bupaF, a, a, a, adag, bupc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-sqrt(8) * g, {bupaF, adag, adag, adag, a, bupc}, op_sites, inst_ops, inst_ops_idxs_set);

          //4
          mpo_gen.AddTerm(-g, {bdna, x, m + sqrt(3) * n, m, m, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g, {bdna, x, sqrt(5) * m + sqrt(7) * n, n, m, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, x, sqrt(9) * m + sqrt(11) * n, m, n, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, x, sqrt(13) * m + sqrt(15) * n, n, n, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bdna, a, adag, sqrt(2) * m + sqrt(6) * n, m, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, adag, a, sqrt(2) * m + sqrt(6) * n, m, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, a, adag, sqrt(10) * m + sqrt(14) * n, n, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, adag, a, sqrt(10) * m + sqrt(14) * n, n, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-g,
                          {bdna, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);
          mpo_gen.AddTerm(-g,
                          {bdna, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1, Fbdnc},
                          op_sites,
                          inst_ops,
                          inst_ops_idxs_set);

          mpo_gen.AddTerm(-sqrt(8) * g, {bdna, a, a, a, adag, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          mpo_gen.AddTerm(-sqrt(8) * g, {bdna, adag, adag, adag, a, Fbdnc}, op_sites, inst_ops, inst_ops_idxs_set);
          break;
        }
          break;
        default:cout << "Do not support for Np>4" << endl;
          exit(0);
      }
    }
  }

  bool Perturbation = params.Perturbation;
  double PerturbationAmplitude = params.PA;
  size_t ChargePeriod = params.PerturbationPeriod;
  if (Perturbation) {
    for (size_t i = 0; i < ElectronSite.size(); i++) {
      size_t x = i / Ly;
      double amplitude =
          -PerturbationAmplitude * cos(M_PI / (double) ChargePeriod + (double) x * (2 * M_PI / (double) ChargePeriod));
      mpo_gen.AddTerm(amplitude, nf, ElectronSite[i]);
    }
    cout << "Add perturbation, mu = " << PerturbationAmplitude << endl;
    cout << "Period of perturbation = " << ChargePeriod << endl;
  }

  auto mpo = mpo_gen.Gen();
  cout << "MPO generated." << endl;

  if (!IsPathExist(kMpoPath)) {
    CreatPath(kMpoPath);
  }

  for (size_t i = 0; i < mpo.size(); i++) {
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
    mpo.DumpTen(i, filename);
  }

  qlmps::FiniteMPO<TenElemT, U1U1QN> finite_mpo(mpo);
  finite_mpo.Truncate(1e-16, 1, 1000);

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  return 0;

}