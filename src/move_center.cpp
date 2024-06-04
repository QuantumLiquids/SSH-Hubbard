/**
 * Connonicalize
 * 
 */

#include "gqdouble.h"
#include "qlmps/qlmps.h"
#include <iostream>
#include <fstream>
#include <vector>

#include "myutil.h"
#include "params_case.h"

using namespace qlmps;
using namespace qlten;
using namespace std;
using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;



int Parser(const int argc, char *argv[],
            size_t& from,
            size_t& to,
            size_t& thread);


/**
 * @example ./move_center --from=10, --to=2, --thread=24
 * the previous center is 10, and the center will move to 2;
 */
int main(int argc, char *argv[]){
    size_t from(0), to(0), thread(0);
    Parser(argc, argv,from, to, thread);

    std::cout << "Argument read: "<< std::endl;
    std::cout << "from = " << from << std::endl;
    std::cout << "to = " << to << std::endl;
    std::cout << "thread = " << thread << std::endl;


    qlten::hp_numeric::SetTensorTransposeNumThreads(thread);
    qlten::hp_numeric::SetTensorManipulationThreads(thread);


    const size_t N = GetNumofMps();
    using TenT = Tensor;
    const string temp_path = kRuntimeTempPath;
    const size_t Ly(4), Lx(N/Ly), Np(3);
    std::cout << "Note we set Ly = 4 and Np = 3 by default" << std::endl;

    vector<IndexT>  pb_out_set(N);
    // translation along x(for electron) and translation along y(for electron);
    for(int i =0;i < N; ++i){
        int residue=i%((2*Np+1)*Ly );
        if(residue<(Np+1)*Ly && residue%(Np+1)==0){
            pb_out_set[i] = pb_outF;
        }
        else pb_out_set[i] = pb_outB;
    }

    SiteVec<TenElemT, U1U1QN> sites=SiteVec<TenElemT, U1U1QN>(pb_out_set);

    FiniteMPS<TenElemT, U1U1QN> mps(sites);

    std::string mps_path = kMpsPath;
    mps.LoadTen(from, GenMPSTenName(mps_path, from));
    if(to < from ){
      for(size_t i = from -1 ; i>= to; i--){
        mps.LoadTen(i, GenMPSTenName(mps_path, i));
        mps.RightCanonicalizeTen(i+1);
        std::cout << "Right canonicalized tensor " << i+1  << std::endl;
        mps.DumpTen(
          i+1,
          GenMPSTenName(mps_path, i+1),
          true
        );
      }
    }else if(to > from){
      for(size_t i = from + 1 ; i<= to; i++){
        mps.LoadTen(i, GenMPSTenName(mps_path, i));
        mps.LeftCanonicalizeTen(i-1);
        std::cout << "Left canonicalized tensor " << i-1  << std::endl;
        mps.DumpTen(
          i-1,
          GenMPSTenName(mps_path, i-1),
          true
        );
      }
    }
    mps.DumpTen(
          to,
          GenMPSTenName(mps_path, to),
          true
      );
    return 0;
}

int Parser(const int argc, char *argv[],
            size_t& from,
            size_t& to,
            size_t& thread){
int nOptionIndex = 1;

string arguement1 = "--from=";
string arguement2 = "--to=";
string arguement3 = "--thread=";
bool from_argument_has(false), to_argument_has(false), thread_argument_has(false);
while (nOptionIndex < argc){
  if (strncmp(argv[nOptionIndex], arguement1.c_str() , arguement1.size()) == 0){
    std::string para_string = &argv[nOptionIndex][arguement1.size()];  
    from = atoi(para_string.c_str());
    from_argument_has = true;
  }else if (strncmp(argv[nOptionIndex], arguement2.c_str(), arguement2.size()) == 0){
    std::string para_string = &argv[nOptionIndex][arguement2.size()]; 
    to = atoi(para_string.c_str());
    to_argument_has = true;
  }else if (strncmp(argv[nOptionIndex], arguement3.c_str(), arguement3.size()) == 0 ){
    std::string para_string = &argv[nOptionIndex][arguement3.size()];
    thread = atoi(para_string.c_str());
    thread_argument_has = true;
  }
  else{
    cout << "Options '" << argv[nOptionIndex] << "' not valid. Run '" << argv[0] << "' for details." << endl;
  //   return -1;
  }
  nOptionIndex++;
}

if(!from_argument_has){
  from = 0;
  std::cout << "Note: no from argument, set it as 0 by default."  << std::endl;
}

if(!to_argument_has){
    to = 10;
    std::cout << "Note: no to argument, set it as 10 by default."  << std::endl;
}

if(!thread_argument_has){
    thread = 24;
    std::cout << "Note: no thread argument, set it as 24 by default." << std::endl;
}

return 0;

}