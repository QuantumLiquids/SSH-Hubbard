// Generate the MPO of Hubbard model with EPC on plaquette.
#include "gqdouble.h"
#include "operators.h"
#include <time.h>
#include <vector>
#include <stdlib.h>     // system
#include "gqmps2/gqmps2.h"

using namespace gqmps2;
using namespace gqten;
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

  if(Ly != 4) {
    std::cout << "Unexpected: Ly = " << Ly << ". The program hopes Ly = 4. " << std::endl;
    exit(1);
  }
  clock_t startTime,endTime;
  startTime = clock();

  OperatorInitial();

  vector<IndexT2>  pb_out_set(N);
  vector<long> Tx(N,-1), Ty(N,-1), ElectronSite(Lx*Ly);
  // translation along x(for electron) and translation along y(for electron);
  auto iter = ElectronSite.begin();
  const size_t total_Ly = Ly+Ly/2*Np;
  for(size_t i = 0;i < N; ++i){
    size_t residue=i%total_Ly;
    if(residue == 0 || residue == Np+1 || residue == 2*(Np+1) || residue == 2*(Np+1)+1 ){
      pb_out_set[i] = pb_outF;
      *iter = i;
      iter++;
    }
    else pb_out_set[i] = pb_outB;
  }

  SiteVec<TenElemT, U1U1QN> sites=SiteVec<TenElemT, U1U1QN>(pb_out_set);
  gqmps2::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  for(size_t i=0;i<N;++i) {
    size_t residue=i % total_Ly;
    if(residue == 0 || residue == Np+1 || residue == 2*(Np+1) || residue == 2*(Np+1)+1 ){
      mpo_gen.AddTerm(U, Uterm, i);
      std::cout << "add site" << i << "Hubbard U term" << endl;
      if(i<N-total_Ly) Tx[i] = i+total_Ly;
      if(residue < 2*(Np+1) ) {
        Ty[i] = i+Np+1;
      } else if (residue== 2*(Np+1)){
        Ty[i] = i + 1;
      } else {
        Ty[i] = i - total_Ly + 1;
      }
    } else {
      unsigned residue2 = residue % (Np+1);
      mpo_gen.AddTerm(omega, (TenElemT)(pow(2,residue2-1))*n_a, i);
      cout << "add site" << i << "phonon potential term (phonon number = "
           <<pow(2,residue2-1) <<")" <<endl;
    }
  }

  for (unsigned i = 0; i < N; i++){
    if(Tx[i]!=-1){
      size_t site1(i), site2(Tx[i]);
      vector<size_t> inst_idxs(Ly-1);
      auto iter = find(ElectronSite.cbegin(), ElectronSite.cend(), site1);
      for(size_t j = 1; j<Ly; j++){
        inst_idxs[j-1] = *(iter+j);
      }
      std::cout << "add site (" << site1 << "," << site2 << ")  hopping  term" << endl;
      mpo_gen.AddTerm( -t, bupcF, site1, bupa, site2, f, inst_idxs);
      mpo_gen.AddTerm( -t, bdnc,  site1, Fbdna,site2, f, inst_idxs);
      mpo_gen.AddTerm(  t, bupaF, site1, bupc, site2, f, inst_idxs);
      mpo_gen.AddTerm(  t, bdna,  site1, Fbdnc,site2, f, inst_idxs);
    }
    if(Ty[i]!=-1){
      if(Ty[i] > i){
        unsigned site1(i), site2(Ty[i]);
        cout << "add site (" << site1 << "," << site2 << ")  hopping  term" << endl;
        mpo_gen.AddTerm( -t, bupcF, site1, bupa, site2);
        mpo_gen.AddTerm( -t, bdnc,  site1, Fbdna,site2);
        mpo_gen.AddTerm(  t, bupaF, site1, bupc, site2);
        mpo_gen.AddTerm(  t, bdna,  site1, Fbdnc,site2);
      }else{
        size_t site1(Ty[i]), site2(i);
        vector<size_t> inst_idxs(Ly-2);
        for(size_t j = 1; j<Ly-1; j++) {
          inst_idxs[j-1] = site1 + j*(Np+1);
        }
        std::cout << "add site (" << site1 << "," << site2 << ")  hopping  term" << endl;
        mpo_gen.AddTerm( -t, bupcF, site1, bupa, site2, f, inst_idxs);
        mpo_gen.AddTerm( -t, bdnc,  site1, Fbdna,site2, f, inst_idxs);
        mpo_gen.AddTerm(  t, bupaF, site1, bupc, site2, f, inst_idxs);
        mpo_gen.AddTerm(  t, bdna,  site1, Fbdnc,site2, f, inst_idxs);
      }
    }
  }

  for (size_t i = 0; i < N; i++){
    //electron phonon interactions
    size_t residue=i % total_Ly;
    if(residue == 0 || residue == (Np+1)){
      vector<size_t> op_sites(Np+3);
      for(size_t j = 0; j<Np+2; j++){
        op_sites[j] = i+j;
      }
      op_sites[Np+2] = Ty[Ty[i]];

      switch(Np){
        case 1:
          mpo_gen.AddTerm( g, {bupcF, x, f, bupa }, op_sites);
          mpo_gen.AddTerm( g, {bdnc,  x, f, Fbdna}, op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x, f, bupc }, op_sites);
          mpo_gen.AddTerm(-g, {bdna,  x, f, Fbdnc}, op_sites);
          cout << "sites ( "<<i << "," << i+1 <<"," << i+2 <<") electron-phonon interaction term" <<endl;
          break;
        case 2:{
          auto op1 = (sqrt(3)-1)*n_a+idB+sqrt(2)*a;
          auto op2 = sqrt(2)*(adag+(-a));
          mpo_gen.AddTerm( g, {bupcF, x,  op1, f, bupa },op_sites);
          mpo_gen.AddTerm( g, {bupcF, a,  op2, f, bupa },op_sites);
          mpo_gen.AddTerm( g, {bdnc,  x,  op1, f, Fbdna},op_sites);
          mpo_gen.AddTerm( g, {bdnc,  a,  op2, f, Fbdna},op_sites);
          mpo_gen.AddTerm(-g, {bupaF, x,  op1, f, bupc },op_sites);
          mpo_gen.AddTerm(-g, {bupaF, a,  op2, f, bupc },op_sites);
          mpo_gen.AddTerm(-g, {bdna,  x,  op1, f, Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,  a,  op2, f, Fbdnc},op_sites);
        }break;
        case 3:{
          auto op1 = sqrt(2)*aadag+sqrt(6)*n_a;
          auto op2 = sqrt(3)*aadag+sqrt(7)*n_a;
          auto op3 =  aadag+sqrt(5)*n_a;
          mpo_gen.AddTerm(g,{bupcF,a,   a,2.0*adag, f, bupa},op_sites);
          mpo_gen.AddTerm(g,{bupcF,adag,adag,2.0*a, f, bupa},op_sites);
          mpo_gen.AddTerm(g,{bupcF,a,   adag, op1,  f, bupa},op_sites);
          mpo_gen.AddTerm(g,{bupcF,adag,a,    op1,  f, bupa},op_sites);
          mpo_gen.AddTerm(g,{bupcF,x,   n_a,  op2,  f, bupa},op_sites);
          mpo_gen.AddTerm(g,{bupcF,x,   aadag,op3,  f, bupa},op_sites);

          mpo_gen.AddTerm(g,{bdnc,a,   a,2.0*adag, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g,{bdnc,adag,adag,2.0*a, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g,{bdnc,a,   adag, op1,  f, Fbdna},op_sites);
          mpo_gen.AddTerm(g,{bdnc,adag,a,    op1,  f, Fbdna},op_sites);
          mpo_gen.AddTerm(g,{bdnc,x,   n_a,  op2,  f, Fbdna},op_sites);
          mpo_gen.AddTerm(g,{bdnc,x,   aadag,op3,  f, Fbdna},op_sites);


          mpo_gen.AddTerm(-g,{bupaF,a,   a,2.0*adag, f, bupc},op_sites);
          mpo_gen.AddTerm(-g,{bupaF,adag,adag,2.0*a, f, bupc},op_sites);
          mpo_gen.AddTerm(-g,{bupaF,a,   adag, op1,  f, bupc},op_sites);
          mpo_gen.AddTerm(-g,{bupaF,adag,a,    op1,  f, bupc},op_sites);
          mpo_gen.AddTerm(-g,{bupaF,x,   n_a,  op2,  f, bupc},op_sites);
          mpo_gen.AddTerm(-g,{bupaF,x,   aadag,op3,  f, bupc},op_sites);

          mpo_gen.AddTerm(-g,{bdna,a,   a,2.0*adag, f, Fbdnc},op_sites);
          mpo_gen.AddTerm(-g,{bdna,adag,adag,2.0*a, f, Fbdnc},op_sites);
          mpo_gen.AddTerm(-g,{bdna,a,   adag, op1,  f, Fbdnc},op_sites);
          mpo_gen.AddTerm(-g,{bdna,adag,a,    op1,  f, Fbdnc},op_sites);
          mpo_gen.AddTerm(-g,{bdna,x,   n_a,  op2,  f, Fbdnc},op_sites);
          mpo_gen.AddTerm(-g,{bdna,x,   aadag,op3,  f, Fbdnc},op_sites);
        }break;
        case 4:{
          auto m = aadag;
          auto n = n_a;   auto& P0=m;auto& P1=n;
          //1
          mpo_gen.AddTerm(g, {bupcF,x,m+sqrt(3)*n,          m, m, f, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF,x,sqrt(5)*m+sqrt(7)*n,  n, m, f, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF,x,sqrt(9)*m+sqrt(11)*n, m, n, f, bupa}, op_sites);
          mpo_gen.AddTerm(g, {bupcF,x,sqrt(13)*m+sqrt(15)*n,n, n, f, bupa}, op_sites);

          mpo_gen.AddTerm(g, {bupcF,a,adag,sqrt(2)*m+sqrt(6)*n,  m, f, bupa},op_sites);
          mpo_gen.AddTerm(g, {bupcF,adag,a,sqrt(2)*m+sqrt(6)*n,  m, f, bupa},op_sites);
          mpo_gen.AddTerm(g, {bupcF,a,adag,sqrt(10)*m+sqrt(14)*n,n, f, bupa},op_sites);
          mpo_gen.AddTerm(g, {bupcF,adag,a,sqrt(10)*m+sqrt(14)*n,n, f, bupa},op_sites);

          mpo_gen.AddTerm(g, {bupcF,a,   a,   adag, sqrt(4)*P0+sqrt(12)*P1, f, bupa},op_sites);
          mpo_gen.AddTerm(g, {bupcF,adag,adag,a,    sqrt(4)*P0+sqrt(12)*P1, f, bupa},op_sites);

          mpo_gen.AddTerm(sqrt(8)*g, {bupcF,a,   a,   a,   adag, f, bupa},op_sites);
          mpo_gen.AddTerm(sqrt(8)*g, {bupcF,adag,adag,adag,a,    f, bupa},op_sites);

          //2
          mpo_gen.AddTerm(g, {bdnc,x,m+sqrt(3)*n,          m, m, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g, {bdnc,x,sqrt(5)*m+sqrt(7)*n,  n, m, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g, {bdnc,x,sqrt(9)*m+sqrt(11)*n, m, n, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g, {bdnc,x,sqrt(13)*m+sqrt(15)*n,n, n, f, Fbdna},op_sites);

          mpo_gen.AddTerm(g, {bdnc,a,adag,sqrt(2)*m+sqrt(6)*n,  m, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g, {bdnc,adag,a,sqrt(2)*m+sqrt(6)*n,  m, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g, {bdnc,a,adag,sqrt(10)*m+sqrt(14)*n,n, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g, {bdnc,adag,a,sqrt(10)*m+sqrt(14)*n,n, f, Fbdna},op_sites);

          mpo_gen.AddTerm(g, {bdnc,a,   a,   adag,sqrt(4)*P0+sqrt(12)*P1, f, Fbdna},op_sites);
          mpo_gen.AddTerm(g, {bdnc,adag,adag,a,   sqrt(4)*P0+sqrt(12)*P1, f, Fbdna},op_sites);

          mpo_gen.AddTerm(sqrt(8)*g, {bdnc,a,   a,   a,   adag, f, Fbdna},op_sites);
          mpo_gen.AddTerm(sqrt(8)*g, {bdnc,adag,adag,adag,a,    f, Fbdna},op_sites);

          //3
          mpo_gen.AddTerm(-g, {bupaF,x,m+sqrt(3)*n,m,m, f, bupc},op_sites);
          mpo_gen.AddTerm(-g, {bupaF,x,sqrt(5)*m+sqrt(7)*n,n,m,f,bupc},op_sites);
          mpo_gen.AddTerm(-g, {bupaF,x,sqrt(9)*m+sqrt(11)*n,m,n,f,bupc},op_sites);
          mpo_gen.AddTerm(-g, {bupaF,x,sqrt(13)*m+sqrt(15)*n,n,n,f,bupc},op_sites);

          mpo_gen.AddTerm(-g, {bupaF,a,adag,sqrt(2)*m+sqrt(6)*n,m,f,bupc},op_sites);
          mpo_gen.AddTerm(-g, {bupaF,adag,a,sqrt(2)*m+sqrt(6)*n,m,f,bupc},op_sites);
          mpo_gen.AddTerm(-g, {bupaF,a,adag,sqrt(10)*m+sqrt(14)*n,n,f,bupc},op_sites);
          mpo_gen.AddTerm(-g, {bupaF,adag,a,sqrt(10)*m+sqrt(14)*n,n,f,bupc},op_sites);

          mpo_gen.AddTerm(-g,{bupaF,a,a,adag,sqrt(4)*P0+sqrt(12)*P1,f,bupc},op_sites);
          mpo_gen.AddTerm(-g,{bupaF,adag,adag,a,sqrt(4)*P0+sqrt(12)*P1,f,bupc},op_sites);

          mpo_gen.AddTerm(-sqrt(8)*g, {bupaF,a,a,a,adag,f,bupc},op_sites);
          mpo_gen.AddTerm(-sqrt(8)*g, {bupaF,adag,adag,adag,a,f,bupc},op_sites);

          //4
          mpo_gen.AddTerm(-g, {bdna,x,m+sqrt(3)*n,m,m,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,x,sqrt(5)*m+sqrt(7)*n,n,m,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,x,sqrt(9)*m+sqrt(11)*n,m,n,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,x,sqrt(13)*m+sqrt(15)*n,n,n,f,Fbdnc},op_sites);

          mpo_gen.AddTerm(-g, {bdna,a,adag,sqrt(2)*m+sqrt(6)*n,m,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,adag,a,sqrt(2)*m+sqrt(6)*n,m,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,a,adag,sqrt(10)*m+sqrt(14)*n,n,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,adag,a,sqrt(10)*m+sqrt(14)*n,n,f,Fbdnc},op_sites);

          mpo_gen.AddTerm(-g, {bdna,a,a,adag,sqrt(4)*P0+sqrt(12)*P1,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-g, {bdna,adag,adag,a,sqrt(4)*P0+sqrt(12)*P1,f,Fbdnc},op_sites);

          mpo_gen.AddTerm(-sqrt(8)*g, {bdna,a,a,a,adag,f,Fbdnc},op_sites);
          mpo_gen.AddTerm(-sqrt(8)*g, {bdna,adag,adag,adag,a,f,Fbdnc},op_sites);
        }break;
        default:
          cout << "This progress does not support for Np > 4 cases" << endl;
          exit(0);
      }
    }
  }

  bool Perturbation=params.Perturbation;
  float PerturbationAmplitude=params.PA;
  int ChargePeriod=4;
  if(Perturbation) {
    for(size_t i=0;i<ElectronSite.size();i++) {
      int x = i/ Ly;
      double amplitude = -PerturbationAmplitude*cos(M_PI/ChargePeriod+x*(2*M_PI/ChargePeriod));
      mpo_gen.AddTerm(amplitude,nf,ElectronSite[i]);
    }
    cout<<"Add perturbation, mu = "<< PerturbationAmplitude<<endl;
    cout << "Period of perturbation = " << ChargePeriod << endl;
  }

  auto mpo = mpo_gen.Gen();
  cout << "MPO generated." << endl;

  const std::string kMpoPath = "mpo";
  const std::string kMpoTenBaseName = "mpo_ten";


  if (!IsPathExist(kMpoPath)) {
    CreatPath(kMpoPath);
  }

  for(size_t i=0; i<mpo.size();i++){
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i) + "." + kGQTenFileSuffix;
    mpo.DumpTen(i,filename);
  }


  endTime = clock();
  cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  return 0;


}