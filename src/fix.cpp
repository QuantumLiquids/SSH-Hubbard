/**
 * Connonicalize
 * 
 */

#include "gqdouble.h"
#include "gqmps2/gqmps2.h"
#include <iostream>
#include <fstream>
#include <vector>

#include "myutil.h"
#include "params_case.h"

using namespace gqmps2;
using namespace gqten;
using namespace std;
using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1U1QN>;

int main(int argc, char *argv[]){
    size_t TotalThreads =56;
    gqten::hp_numeric::SetTensorTransposeNumThreads(TotalThreads);
    gqten::hp_numeric::SetTensorManipulationTotalThreads(TotalThreads);
    

    CaseParams params(argv[1]);
    
    size_t Lx = params.Lx, Ly = params.Ly, Np = params.Np;
    size_t N = Lx * Ly + (2 * Lx * Ly - Ly) * Np;
    if(GetNumofMps()!=N){
        std::cout << "The number of mps files are inconsistent with mps size!" <<std::endl;
        exit(1);
    }

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


    SiteVec<TenElemT, U1U1QN> sites=SiteVec<TenElemT, U1U1QN>(pb_out_set);
    FiniteMPST mps(sites);
    mps.Load();
    cout << "mps loaded" <<endl;
    mps.Centralize(0);
    mps.Dump();

    return 0;


}