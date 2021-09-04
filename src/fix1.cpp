/**
 * Fix single site update, svdinfo==0 collepse
 * 
 */

#include "gqdouble.h"
// #include "gqmps2/algorithm/lanczos_solver.h"
// #include "gqmps2/algorithm/lanczos_solver_impl.h"
#include "gqmps2/gqmps2.h"
#include <iostream>
#include <fstream>
#include <vector>
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using namespace gqmps2;
using namespace gqten;
using namespace std;

int main(){
    std::cout << "This program used to patch a left site mps by single site Lanczos" <<std::endl;
    std::cout << "The input must include these files: renv.gqten, lenv.gqten, mpo_ten_r.gqten, mpo_ten_l.gqten, mps_ten_r.gqten" <<std::endl;
    std::cout << "The output is the file mps_ten_l.gqten" <<std::endl;
    std::cout << "where lenv.gqten is mps_ten_l.gqten's left environment" <<std::endl;
    std::cout << "and renv.gqten is mps_ten_r.gqten's right environment" <<std::endl;
    size_t TotalThreads =56;
    gqten::hp_numeric::SetTensorTransposeNumThreads(TotalThreads);
    gqten::hp_numeric::SetTensorManipulationTotalThreads(TotalThreads);
    
    size_t N = 1332; //total site number
    size_t target_site = 567; 
    size_t lenv = target_site;
    size_t renv = (N-1) - target_site;

    Tensor renv_ten_before, renv_ten, lenv_ten, mpo1, mpo2, mps1, mps2;
    //mps1 is the target tensor
    string  file = "renv.gqten";
    ifstream tensor_file(file, ifstream::binary);
    tensor_file >> renv_ten_before;
    tensor_file.close();

    file = "lenv.gqten";
    tensor_file.open(file, ifstream::binary);
    tensor_file >> lenv_ten;
    tensor_file.close();

    file = "mps_ten_r.gqten";
    tensor_file.open(file, ifstream::binary);
    tensor_file >> mps2;
    tensor_file.close();

    file = "mpo_ten_r.gqten";
    tensor_file.open(file, ifstream::binary);
    tensor_file >> mpo2;
    tensor_file.close();

    file = "mpo_ten_l.gqten";
    tensor_file.open(file, ifstream::binary);
    tensor_file >> mpo1;
    tensor_file.close();

    Tensor mps2_dag = Dag(mps2);
    Tensor temp1, temp2;
    Contract(&mps2, &renv_ten_before, {{2},{0}}, &temp1);
    Contract(&temp1, &mpo2, {{1,2},{1,3}},&temp2);
    Contract(&temp2, &mps2_dag, {{3,1},{1,2}}, &renv_ten);
    temp1 = Tensor();
    temp2 = Tensor();

    IndexT2 index0 = InverseIndex( lenv_ten.GetIndexes()[0] ) ;
    IndexT2 index1 = InverseIndex( mpo1.GetIndexes()[1] ) ;
    IndexT2 index2 = InverseIndex( renv_ten.GetIndexes()[0] ) ;
    vector<IndexT2> indexes = {index0, index1, index2};
    Tensor* initial_state = new Tensor({index0, index1, index2});

    gqten::ShapeT blk_shape = {index0.GetQNSctNum(),
                        index1.GetQNSctNum(),
                        index2.GetQNSctNum() };
    gqten::CoorsT blk_coors;
    bool flag(false);
    for(size_t i =0; i< blk_shape[0];i++){
        for(size_t j=0;j<blk_shape[1];j++){
            for(size_t k =0;k<blk_shape[2];j++){
                if (CalcDiv(indexes, {i,j,k }) == qn0) {
                    blk_coors = {i,j,k};
                    flag = true;
                    break;
                }
            }
            if(flag){
                break;
            }
        }
        if(flag){
            break;
        }
    }
    gqten::CoorsT zeros_coor = {0,0,0};
    gqten::BlockSparseDataTensor<gqten::GQTEN_Double, U1U1QN>& bstd = initial_state->GetBlkSparDataTen();
    bstd.ElemSet(std::make_pair(blk_coors,zeros_coor), 1.0);

    gqmps2::LanczosParams params(1e-7, 200);

    std::vector<Tensor *> eff_ham(3);
    eff_ham[0] = const_cast<Tensor *>(&lenv_ten);
    eff_ham[1] = const_cast<Tensor *>(&mpo1);
    eff_ham[2] = const_cast<Tensor *>(&renv_ten);

    LanczosRes<Tensor> res = LanczosSolver<Tensor>(
                            eff_ham, 
                            initial_state,
                            &eff_ham_mul_single_site_state, 
                            params);
    mps1 = std::move(*res.gs_vec);
    std::cout << "ground state energy = " << res.gs_eng <<std::endl;
    delete res.gs_vec;
    
    file = "mps_ten_l.gqten";
    ofstream dump_file(file, ofstream::binary);
    dump_file << mps1;

    return 0;


}