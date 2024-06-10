//DMRG variation with only one process
#include "gqdouble.h"
#include "operators.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "qlmps/qlmps.h"
#include "singlesiteupdate2.h"
#include "twositeupdate2.h"
#include "myutil.h"


using namespace qlmps;
using namespace qlten;
using namespace std;


#include "params_case.h"


int main(int argc, char *argv[]) {
    CaseParams params(argv[1]);

    unsigned Lx=params.Lx, Ly = params.Ly, Np = params.Np;
    unsigned N = (4+2*Np)*Lx;
    cout << "System size = (" << Lx << "," <<Ly<<")"<< endl;
    cout << "The number of electron sites =" <<Lx*Ly<<endl;
    cout << "The number of phonon pseudosite (per bond) =" << Np <<endl;
    cout << "The number of phonon pseudosite (total) =" << (2*Lx*Ly-Ly)*Np <<endl;
    cout << "The total number of sites = " << N<<endl;
    float t = params.t, g = params.g, U = params.U, omega = params.omega;
    cout << "Model parameter: t =" << t <<", g =" << g <<", U =" << U <<",omega=" << omega<< endl;

    clock_t startTime,endTime;
    startTime = clock();
    if(Ly != 4) {
      std::cout << "Unexpected: Ly = " << Ly << ". The program hopes Ly = 4. " << std::endl;
      exit(1);
    }

    OperatorInitial();

    vector<IndexT2>  pb_out_set(N);
    vector<long> Tx(N,-1), Ty(N,-1), ElectronSite(Lx*Ly);
    auto iter = ElectronSite.begin();
    const size_t total_Ly = Ly+Ly/2*Np;
    for(size_t i =0;i < N; ++i) {
      size_t residue=i%total_Ly;
      if(residue == 0 || residue == Np+1 || residue == 2*(Np+1) || residue == 2*(Np+1)+1 ){
        pb_out_set[i] = pb_outF;
        *iter = i;
        iter++;
      }
      else pb_out_set[i] = pb_outB;
    }

    SiteVec<TenElemT, U1U1QN> sites=SiteVec<TenElemT, U1U1QN>(pb_out_set);

    MPO<Tensor> mpo(N);
    const std::string kMpoPath = "mpo";
    const std::string kMpoTenBaseName = "mpo_ten";
    for(size_t i=0; i<mpo.size();i++){
        std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
        mpo.LoadTen(i,filename);
    }
   
    cout << "MPO loaded." << endl;
    using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
    FiniteMPST mps(sites);

    std::vector<long unsigned int> stat_labs(N,1);
    int sitenumber_perhole;
    if(params.Numhole>0){
        sitenumber_perhole = ElectronSite.size()/params.Numhole;
    }else{
        sitenumber_perhole = 4*ElectronSite.size();
    }

    srand((unsigned)time(NULL));
    int qn_label = 1;
    for (int i = 0; i <ElectronSite.size() ; i++) {
        if(i%sitenumber_perhole==sitenumber_perhole/2){
            stat_labs[ElectronSite[i]] = 3;
        }else{
            stat_labs[ElectronSite[i]] = qn_label;
            qn_label=3-qn_label;
        }
    }


    qlten::hp_numeric::SetTensorManipulationThreads(params.TotalThreads);
    qlmps::SingleVMPSSweepParams sweep_params(
        params.Sweeps,
        params.Dmin, params.Dmax, params.CutOff,
        qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
        params.noise
    );
    if(IsPathExist(kMpsPath)){
        if(N== GetNumofMps() ){
            cout << "The number of mps files is consistent with mps size." <<endl;
            cout << "Directly use mps from files." <<endl;
        }else{
            qlmps::DirectStateInitMps(mps, stat_labs);
            cout << "Initial mps as direct product state." <<endl;
            mps.Dump(sweep_params.mps_path, true);
        }
    }else{
        qlmps::DirectStateInitMps(mps, stat_labs);
        cout << "Initial mps as direct product state." <<endl;
        mps.Dump(sweep_params.mps_path, true);
    }
    auto e0 = qlmps::TwoSiteFiniteVMPS2(mps, mpo, sweep_params);
    std::cout << "E0/site: " << e0 / N << std::endl;

    endTime = clock();
    cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    return 0;


}
