#include "gqmps2/gqmps2.h"
using gqmps2::kMpsPath;
using gqmps2::kMpsTenBaseName;
using gqmps2::kGQTenFileSuffix;

//number of mps file in default mps path("./mps")
size_t GetNumofMps(){
    size_t NumberOfMpsFile = 0;
    for(NumberOfMpsFile=0; NumberOfMpsFile < 1e5; NumberOfMpsFile++ ){
        std::string file;
        file = kMpsPath + "/" + kMpsTenBaseName + std::to_string(NumberOfMpsFile) + "." + kGQTenFileSuffix;
        std::ifstream ifs(file, std::ifstream::binary);
        if(ifs.good()){
            ifs.close();
        }else{
            break;
        }
    }
    return NumberOfMpsFile;
}


//If site `i` is electron site, for SSH-Hubbard model
bool IsElectron(size_t i, size_t Ly, size_t Np){
    int residue=i%((2*Np+1)*Ly );
    if(residue<(Np+1)*Ly && residue%(Np+1)==0){
        return true;
    }
    else return false;
}


void Show(std::vector<size_t> v){
    for(auto iter = v.begin(); iter<v.end(); iter++){
        std::cout<< *iter<<",";
    }
    std::cout <<'\b' <<std::endl;
}
