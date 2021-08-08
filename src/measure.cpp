/*
    measure.cpp
    for measure one point and two point function
*/
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"

#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include <time.h>
#include <stdlib.h>

#include "myutil.h"
#include "my_measure.h"

#include "gqten/utility/timer.h"  

using std::cout;
using std::endl;
using std::vector;
using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1U1QN>;
using gqmps2::SiteVec;
using gqmps2::MeasureOneSiteOp;
using gqten::Timer;

int main(int argc, char *argv[]) {
    clock_t startTime, endTime;
    startTime = clock();


    CaseParams params(argv[1]);
    
    size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
    size_t N = Lx * Ly + (2 * Lx * Ly - Ly) * Np;
    if(GetNumofMps()!=N){
        std::cout << "The number of mps files are inconsistent with mps size!" <<std::endl;
        exit(1);
    }

    OperatorInitial();

    std::vector<IndexT2>  pb_out_set(N);
    std::vector<size_t> Fsite_set(Lx*Ly);
    std::vector<size_t> Bsite_set(N-Lx*Ly);
    auto iterF = Fsite_set.begin();
    auto iterB = Bsite_set.begin();

     for(size_t i =0;i < N; ++i){
        if(IsElectron(i, Ly, Np)){
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
    mps.Load();
    cout << "mps loaded" <<endl;
    gqten::hp_numeric::SetTensorTransposeNumThreads(params.TotalThreads);
    gqten::hp_numeric::SetTensorManipulationTotalThreads(params.TotalThreads);
    gqten::hp_numeric::SetTensorDecompOuterParallelThreads(params.SvdOuterThreads);
    
    cout << "bond dimension of middel mps = " ;
    cout << mps[N/2].GetShape()[0] <<endl;


    Timer onesite_timer("measure one site operators");
    MeasureOneSiteOp(mps, {sx, sz, nf}, Fsite_set, {"sx","sz","nf"});
    MeasureOneSiteOp(mps, n_a, Bsite_set, "nphonon");
    cout << "measured one point function.<====" <<endl;
    onesite_timer.PrintElapsed();

    vector<vector<size_t>> two_point_sites_setF;
    
    size_t beginx = Lx/4;
//    int beginx = 4; //for Lx=24
    size_t endx = beginx+Lx/2;
    two_point_sites_setF.reserve(Ly*(endx-beginx));
    for (size_t y = 0; y < Ly; ++y) {
        size_t site1F = beginx*Ly+y;
        for (size_t x = beginx + 1; x < endx; ++x) {
            size_t site2F = x*Ly +y;
            two_point_sites_setF.push_back({Fsite_set[site1F],Fsite_set[site2F]});
        }
    }

    Timer twosite_timer("measure two site operators");
    MeasureTwoSiteOp(mps, sz,   sz,  two_point_sites_setF, "szsz");
    MeasureTwoSiteOp(mps, sx,   sx,  two_point_sites_setF, "sxsx");
    MeasureTwoSiteOp(mps, -isy, isy, two_point_sites_setF, "sysy");
    MeasureTwoSiteOp(mps, nf,   nf,  two_point_sites_setF, "nfnf");
    MeasureTwoSiteOp(mps, cupccdnc,cdnacupa,two_point_sites_setF, "onsitesc");
    cout << "measured two point function.<====" <<endl;
    twosite_timer.PrintElapsed();

    endTime = clock();
    cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    return 0;

}

