/*
    measure.cpp
    for measure one point function. memory optimized version
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

int main(int argc, char *argv[]) {
    clock_t startTime, endTime;
    startTime = clock();

    CaseParams params(argv[1]);
    
    size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
    size_t N = Lx*Ly+(2*Lx*Ly - Ly - Lx)*Np;
    if(GetNumofMps()!=N){
        std::cout << "The number of mps files are inconsistent with mps size!" <<std::endl;
        exit(1);
    }

    OperatorInitial();

    std::vector<IndexT>  pb_out_set(N);
    std::vector<size_t> Fsite_set(Lx*Ly);
    std::vector<size_t> Bsite_set(N-Lx*Ly);
    auto iterF = Fsite_set.begin();
    auto iterB = Bsite_set.begin();


    const size_t total_Ly = (2*Np+1)*Ly - Np;
    for(size_t i =0;i < N; ++i){
      size_t residue=i%total_Ly;
      if(residue<(Np+1)*Ly && residue%(Np+1)==0){
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
    qlten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads);
    qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);

    Timer onesite_timer("measure one site operators");
    MeasureOneSiteOp(mps, {sz, nf}, Fsite_set, {"sz","nf"});
    MeasureOneSiteOp(mps, n_a, Bsite_set, "nphonon");
    cout << "measured one point function.<====" <<endl;
    onesite_timer.PrintElapsed();





    endTime = clock();
    cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    return 0;

}

