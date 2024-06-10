/*
    oy_measureSC_pla.cpp
    for measure pair correlation function. For plaquette symmetry. Phi_xx
    Delta_yx
    usage:
        mpirun -n 4 ./ox_measureSC_ps2
    note: processor number must be 4.
    Optional arguments:
      --start=
      --end=
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

#include "boost/mpi.hpp"

using std::cout;
using std::endl;
using std::vector;
using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
using qlmps::SiteVec;
using qlmps::MeasureOneSiteOp;
using qlten::Timer;
using qlmps::MeasureElectronPhonon4PointFunction;


int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;
  clock_t startTime, endTime;
  startTime = clock();

  size_t beginx;
  size_t endx;
  bool start_argument_has = Parser(argc, argv, beginx, endx);

  CaseParams params(argv[1]);

  size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
  size_t N = Lx*Ly+(Lx*Ly-Ly)*Np;
  if(GetNumofMps()!=N){
    std::cout << "The number of mps files are inconsistent with mps size!" <<std::endl;
    exit(1);
  }

  if( !start_argument_has ) {
    beginx = Lx/4;
    endx = beginx+Lx/2;
  }

  OperatorInitial();

  std::vector<IndexT2>  pb_out_set(N);
  std::vector<size_t> Fsite_set(Lx*Ly);
  std::vector<size_t> Bsite_set(N-Lx*Ly);
  auto iterF = Fsite_set.begin();
  auto iterB = Bsite_set.begin();

  const size_t total_Ly = Ly * (Np + 1);
  for(size_t i =0;i < N; ++i){
    size_t residue=i%total_Ly;
    if(residue < Ly) {
      pb_out_set[i] = pb_outF;
      *iterF = i;
      iterF++;
    }
    else {
      pb_out_set[i] = pb_outB;
      *iterB = i;
      iterB++;
    }
  }

  cout << "The Fermion sites: "<<endl;
  Show(Fsite_set);
  cout << '\n';
  cout << "The Boson sites:" <<endl;
  Show(Bsite_set);

  SiteVec<TenElemT, U1U1QN> sites=SiteVec<TenElemT, U1U1QN>(pb_out_set);
  FiniteMPST mps(sites);

  qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);


  Timer foursite_timer("measure pair-pair correlation for justifying the pairing symmetry. ");
  vector<vector<size_t>> fourpoint_sitessetF; // Phi_xx
  std::vector<size_t> Tx(Lx*Ly), Ty(Lx*Ly);
  for (size_t i = 0; i < Lx*Ly; ++i){
    size_t y = i%Ly, x = i/Ly;
    Tx[i] = y +Ly*((x+1)%Lx);
    Ty[i] = (y+1)%Ly + Ly*x;
  }

  fourpoint_sitessetF.reserve(Ly*(endx-beginx));


  for (size_t y = 0; y < Ly; ++y){
    auto site1F = beginx*Ly+y;
    for (size_t x = beginx+2; x < endx; x=x+1) {
      auto site2F = x * Ly + y;
      vector<size_t> xxsites={Fsite_set[site1F], Fsite_set[Tx[site1F]],Fsite_set[site2F],Fsite_set[Tx[site2F]]};
      sort(xxsites.begin(),xxsites.end());
      fourpoint_sitessetF.push_back(xxsites);
    }
  }


  std::vector<Tensor> sc_phys_ops_a = { bupcF, Fbdnc, bupaF, Fbdna };
  std::vector<Tensor> sc_phys_ops_b = { bdnc, bupc, bupaF, Fbdna };
  std::vector<Tensor> sc_phys_ops_c = { bupcF, Fbdnc, bdna, bupa };
  std::vector<Tensor> sc_phys_ops_d = { bdnc, bupc, bdna, bupa };

  std::string file_name_postfix;
  if (start_argument_has) {
    file_name_postfix = "begin" + std::to_string(beginx) + "end" + std::to_string(endx);
  } else {
    file_name_postfix = "";
  }

  if(world.rank()==0){
    MeasureElectronPhonon4PointFunction(mps, sc_phys_ops_a, fourpoint_sitessetF, Np, "scsxxa" + file_name_postfix);//PS for plaquette symmetry
  }else if(world.rank()==1){
    MeasureElectronPhonon4PointFunction(mps, sc_phys_ops_b, fourpoint_sitessetF, Np, "scsxxb" + file_name_postfix);
  }else if(world.rank()==2){
    MeasureElectronPhonon4PointFunction(mps, sc_phys_ops_c, fourpoint_sitessetF, Np, "scsxxc" + file_name_postfix);
  }else if(world.rank()==3){
    MeasureElectronPhonon4PointFunction(mps, sc_phys_ops_d, fourpoint_sitessetF, Np, "scsxxd" + file_name_postfix);
  }

  cout << "measured x--x four point function.<====" <<endl;
  foursite_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  return 0;

}
