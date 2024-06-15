#include "gqdouble.h"
#include "operators.h"
#include <vector>
#include "qlmps/qlmps.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

#include "params_case.h"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <params file> <phonon x file>" << std::endl;
    return 1;
  }
  CaseParams params(argv[1]);

  const size_t Lx = params.Lx, Ly = params.Ly;
  const size_t N = Lx * Ly;
  cout << "System size = (" << Lx << "," << Ly << ")" << endl;
  cout << "The number of electron sites =" << Lx * Ly << endl;
  cout << "The total number of sites = " << N << endl;
  float t = params.t, g = params.g, U = params.U, omega = params.omega;
  cout << "Model parameter: t =" << t << ", g =" << g << ", U =" << U << ",omega=" << omega << endl;

  const char *filename = argv[2];
  std::cout << "phonon x filename : " << filename << std::endl;
  std::ifstream ifs(filename, std::ofstream::binary);
  if (!ifs) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return 1;
  }
  size_t array_size = N + (Lx - 1) * Ly;
  double *x_vec = new double[array_size];
  ifs.read(reinterpret_cast<char *>(x_vec), array_size * sizeof(double));
  if (!ifs) {
    std::cerr << "Error reading file: " << filename << std::endl;
    delete[] x_vec;
    return 1;
  }
  ifs.close();

  clock_t startTime, endTime;
  startTime = clock();


  // average reflection symmetry
  std::vector<double> x_vec_reverse(x_vec_reverse);
  for (size_t i = 0; i < N + (Lx - 1) * Ly; i++) {
    x_vec_reverse[i] = x_vec[i];
  }
  reverse(x_vec_reverse.begin(), x_vec_reverse.end());
  std::cout << "Displacement X : " << std::endl;
  std::cout << "[";
  for (size_t i = 0; i < N + (Lx - 1) * Ly; i++) {
    std::cout << " " << x_vec[i];
  }
  std::cout << "]\n" << std::endl;
  std::cout << "Reflection Averaged Displacement X : ";
  for (size_t i = 0; i < N + (Lx - 1) * Ly; i++) {
    x_vec[i] = (x_vec[i] + x_vec_reverse[i]) / 2;
  }
  std::cout << "[";
  for (size_t i = 0; i < N + (Lx - 1) * Ly; i++) {
    std::cout << " " << x_vec[i];
  }
  std::cout << "]\n" << std::endl;
  // average translation symmetry
  std::vector<double> horizontal_x((Lx - 1)), vertical_x(Lx);
  for (size_t x = 0; x < Lx; x++) {
    size_t start_ph_site = x * 2 * Ly;
    double sum_a = 0.0;
    for (size_t i = start_ph_site; i < start_ph_site + Ly; i++) {
      sum_a += x_vec[i];
    }
    vertical_x[x] = sum_a / Ly;
    if (x < Lx - 1) {
      double sum_b = 0.0;
      for (size_t i = start_ph_site + Ly; i < start_ph_site + 2 * Ly; i++) {
        sum_b += x_vec[i];
      }
      horizontal_x[x] = sum_b / Ly;
    }
  }

  OperatorInitial();
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_outF);
  qlmps::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  for (size_t i = 0; i < N; ++i) {
    mpo_gen.AddTerm(U, Uterm, i);
//    cout << "add site" << i << "Hubbard U term" << endl;
  }

  //horizontal interaction
  for (size_t i = 0; i < N - Ly; ++i) {
    size_t site1 = i, site2 = i + Ly;
    size_t x = i / Ly;
    mpo_gen.AddTerm(-t + g * horizontal_x[x], bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-t + g * horizontal_x[x], bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(t - g * horizontal_x[x], bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(t - g * horizontal_x[x], bdna, site1, Fbdnc, site2, f);
//    cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
  }
  //vertical interaction
  for (size_t i = 0; i < N; ++i) {
    size_t y = i % Ly, x = i / Ly;
    if (y < Ly - 1) {
      size_t site1 = i, site2 = i + 1;
      mpo_gen.AddTerm(-t + g * vertical_x[x], bupcF, site1, bupa, site2);
      mpo_gen.AddTerm(-t + g * vertical_x[x], bdnc, site1, Fbdna, site2);
      mpo_gen.AddTerm(t - g * vertical_x[x], bupaF, site1, bupc, site2);
      mpo_gen.AddTerm(t - g * vertical_x[x], bdna, site1, Fbdnc, site2);
//      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    } else if (Ly > 2) {
      size_t site1 = i - Ly + 1, site2 = i;
      mpo_gen.AddTerm(-t + g * vertical_x[x], bupcF, site1, bupa, site2, f);
      mpo_gen.AddTerm(-t + g * vertical_x[x], bdnc, site1, Fbdna, site2, f);
      mpo_gen.AddTerm(t - g * vertical_x[x], bupaF, site1, bupc, site2, f);
      mpo_gen.AddTerm(t - g * vertical_x[x], bdna, site1, Fbdnc, site2, f);
//      cout << "add site (" << site1 << "," << site2 << ")  hopping term" << endl;
    }
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

  delete[] x_vec;
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  return 0;

}