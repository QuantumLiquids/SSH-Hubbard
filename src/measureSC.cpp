/*
    measureSC.cpp
    for measure pair correlation function
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
using gqmps2::MeasureElectronPhonon4PointFunction;

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
    
    cout << "bond dimension of middel mps = " ;
    cout << mps[N/2].GetShape()[0] <<endl;


    Timer foursite_timer("measure four site operators");
    vector<vector<size_t>> xx_fourpoint_sitessetF;
    vector<vector<size_t>> yy_fourpoint_sitessetF;
    std::vector<size_t> Tx(Lx*Ly), Ty(Lx*Ly);
    for (size_t i = 0; i < Lx*Ly; ++i){
        size_t y = i%Ly, x = i/Ly;
        Tx[i] = y +Ly*((x+1)%Lx);
        Ty[i] = (y+1)%Ly + Ly*x;
    }
    size_t beginx = Lx/4;
//    int beginx = 4; //for Lx=24
    size_t endx = beginx+Lx/2;
    xx_fourpoint_sitessetF.reserve(Ly*(endx-beginx));
    yy_fourpoint_sitessetF.reserve(Ly*(endx-beginx));
    for (size_t y = 0; y < Ly; ++y){
        auto site1F = beginx*Ly+y;
        for (size_t x = beginx+2; x < endx; x=x+1) {
            auto site2F = x * Ly + y;
            vector<size_t> xxsites={Fsite_set[site1F], Fsite_set[Tx[site1F]],Fsite_set[site2F],Fsite_set[Tx[site2F]]};
            xx_fourpoint_sitessetF.push_back(xxsites);


            vector<size_t> yysites={Fsite_set[site1F], Fsite_set[Ty[site1F]],Fsite_set[site2F],Fsite_set[Ty[site2F]]};
	        sort(yysites.begin(),yysites.end());
            yy_fourpoint_sitessetF.push_back(yysites);
        }
    }

    std::vector<Tensor> sc_phys_ops_a = { bupcF, Fbdnc, bupaF, Fbdna };
    std::vector<Tensor> sc_phys_ops_b = { bdnc, bupc, bupaF, Fbdna };
    std::vector<Tensor> sc_phys_ops_c = { bupcF, Fbdnc, bdna, bupa };
    std::vector<Tensor> sc_phys_ops_d = { bdnc, bupc, bdna, bupa };


    MeasureElectronPhonon4PointFunction(mps,sc_phys_ops_a,yy_fourpoint_sitessetF,Np,"scsyya");
    MeasureElectronPhonon4PointFunction(mps,sc_phys_ops_b,yy_fourpoint_sitessetF,Np,"scsyyb");
    MeasureElectronPhonon4PointFunction(mps,sc_phys_ops_c,yy_fourpoint_sitessetF,Np,"scsyyc");
    MeasureElectronPhonon4PointFunction(mps,sc_phys_ops_d,yy_fourpoint_sitessetF,Np,"scsyyd");
    
    cout << "measured y--y four point function.<====" <<endl;
    foursite_timer.PrintElapsed();



    endTime = clock();
    cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    return 0;

}

