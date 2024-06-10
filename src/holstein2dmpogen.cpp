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

  cout << "Generate the MPO for Holstein-Hubbard model" << std::endl;
  size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
  size_t N = (1 + Np) * (Lx * Ly);
  cout << "System size = (" << Lx << "," << Ly << ")" << endl;
  cout << "The number of electron sites =" << Lx * Ly << endl;
  cout << "The number of phonon pseudosite (per bond) =" << Np << endl;
  cout << "The number of phonon pseudosite (total) =" << Np * Lx * Ly << endl;
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
    if (i % (Np + 1) == 0) {
      pb_out_set[i] = pb_outF;
    } else pb_out_set[i] = pb_outB;
  }

  SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(pb_out_set);
  qlmps::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  auto iter = ElectronSite.begin();
  for (size_t i = 0; i < N; ++i) {
    if (i % (Np + 1) == 0) {//fermion site
      size_t fermion_site_num = i / (Np + 1);
      size_t x = fermion_site_num / Ly;
      size_t y = fermion_site_num % Ly;

      mpo_gen.AddTerm(U, Uterm, i);
      std::cout << "add site" << i << "Hubbard U term" << endl;
      *iter = i;
      iter++;
      if (x < Lx - 1) Tx[i] = i + (Np + 1) * Ly;
      if (y != Ly - 1) Ty[i] = i + (Np + 1);
      else Ty[i] = i - (Np + 1) * (Ly - 1);
    } else {//boson site
      size_t residue = i % (Np + 1);
      mpo_gen.AddTerm(omega, (TenElemT) (pow(2, residue - 1)) * n_a, i);
      cout << "add site" << i << "phonon potential term (phonon number = "
           << pow(2, residue - 1) << ")" << endl;
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

  //electron-phonon interaction
  for (size_t i = 0; i < N; i++) {
    if (i % (Np + 1) == 0) {

      std::vector<size_t> op_sites(Np + 1);
      for (size_t j = 0; j < Np + 1; j++) {
        op_sites[j] = i + j;
      }
      switch (Np) {
        case 1:mpo_gen.AddTerm(g, {nf, x}, op_sites);
          cout << "sites ( " << i << "," << i + 1 << ") electron-phonon interaction term" << endl;
          break;
        case 2: {
          auto op1 = (sqrt(3) - 1) * n_a + idB + sqrt(2) * a;
          auto op2 = sqrt(2) * (adag + (-a));
          mpo_gen.AddTerm(g, {nf, x, op1}, op_sites);
          mpo_gen.AddTerm(g, {nf, a, op2}, op_sites);
        }
          break;
        case 3: {
          auto op1 = sqrt(2) * aadag + sqrt(6) * n_a;
          auto op2 = sqrt(3) * aadag + sqrt(7) * n_a;
          auto op3 = aadag + sqrt(5) * n_a;
          mpo_gen.AddTerm(g, {nf, a, a, 2.0 * adag}, op_sites);
          mpo_gen.AddTerm(g, {nf, adag, adag, 2.0 * a}, op_sites);
          mpo_gen.AddTerm(g, {nf, a, adag, op1}, op_sites);
          mpo_gen.AddTerm(g, {nf, adag, a, op1}, op_sites);
          mpo_gen.AddTerm(g, {nf, x, n_a, op2}, op_sites);
          mpo_gen.AddTerm(g, {nf, x, aadag, op3}, op_sites);
        }
          break;
        case 4: {
          auto m = aadag;
          auto n = n_a;
          auto &P0 = m;
          auto &P1 = n;
          mpo_gen.AddTerm(g, {nf, x, m + sqrt(3) * n, m, m}, op_sites);
          mpo_gen.AddTerm(g, {nf, x, sqrt(5) * m + sqrt(7) * n, n, m}, op_sites);
          mpo_gen.AddTerm(g, {nf, x, sqrt(9) * m + sqrt(11) * n, m, n}, op_sites);
          mpo_gen.AddTerm(g, {nf, x, sqrt(13) * m + sqrt(15) * n, n, n}, op_sites);

          mpo_gen.AddTerm(g, {nf, a, adag, sqrt(2) * m + sqrt(6) * n, m}, op_sites);
          mpo_gen.AddTerm(g, {nf, adag, a, sqrt(2) * m + sqrt(6) * n, m}, op_sites);
          mpo_gen.AddTerm(g, {nf, a, adag, sqrt(10) * m + sqrt(14) * n, n}, op_sites);
          mpo_gen.AddTerm(g, {nf, adag, a, sqrt(10) * m + sqrt(14) * n, n}, op_sites);

          mpo_gen.AddTerm(g, {nf, a, a, adag, sqrt(4) * P0 + sqrt(12) * P1}, op_sites);
          mpo_gen.AddTerm(g, {nf, adag, adag, a, sqrt(4) * P0 + sqrt(12) * P1}, op_sites);

          mpo_gen.AddTerm(sqrt(8) * g, {nf, a, a, a, adag}, op_sites);
          mpo_gen.AddTerm(sqrt(8) * g, {nf, adag, adag, adag, a}, op_sites);
        }
          break;
        default:cout << "This progress does not support for Np > 4 cases" << endl;
          exit(0);
      }
    }
  }

  bool Perturbation = params.Perturbation;
  float PerturbationAmplitude = params.PA;
  int ChargePeriod = params.PerturbationPeriod;
  if (Perturbation) {
    for (size_t i = 0; i < ElectronSite.size(); i++) {
      int x = i / Ly;
      double amplitude = -PerturbationAmplitude * cos(M_PI / ChargePeriod + x * (2 * M_PI / ChargePeriod));
      mpo_gen.AddTerm(amplitude, nf, ElectronSite[i]);
    }
    cout << "Add perturbation, mu = " << PerturbationAmplitude << endl;
    cout << "Period of perturbation = " << ChargePeriod << endl;
  }

  auto mpo = mpo_gen.Gen();
  cout << "MPO generated." << endl;

  const std::string kMpoPath = "mpo";
  const std::string kMpoTenBaseName = "mpo_ten";

  if (!IsPathExist(kMpoPath)) {
    CreatPath(kMpoPath);
  }

  for (size_t i = 0; i < mpo.size(); i++) {
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
    mpo.DumpTen(i, filename);
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  return 0;

}