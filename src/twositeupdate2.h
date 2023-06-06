// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
*         Rongyang Sun <sun-rongyang@outlook.com>
* Creation Date: 2021-7-30
*
* Description: GraceQ/MPS2 project. Implementation details for noised two-site algorithm.
*/

/**
@file two_site_update_noise_finite_vmps_impl.h
@brief Implementation details for noised two-site algorithm.
*/

#pragma once


#include "gqmps2/algorithm/vmps/single_site_update_finite_vmps.h"   // SingleVMPSSweepParams
#include "gqmps2/algorithm/vmps/two_site_update_finite_vmps.h"      // helper functions
#include "gqmps2/one_dim_tn/mpo/mpo.h"                              // MPO
#include "gqmps2/one_dim_tn/mps/finite_mps/finite_mps.h"            // FiniteMPS
#include "gqmps2/utilities.h"                                       // IsPathExist, CreatPath
#include "gqmps2/one_dim_tn/framework/ten_vec.h"                    // TenVec
#include "gqmps2/consts.h"
#include "gqten/gqten.h"
#include "gqten/utility/timer.h"                                    // Timer
#include "singlesiteupdate2.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <stdio.h>    // remove
#ifdef Release
#define NDEBUG
#endif
#include <assert.h>

namespace gqmps2{
using namespace gqten;


using TwoSiteVMPSSweepParams = SingleVMPSSweepParams;

// Forward declarition
template <typename DTenT>
inline double MeasureEE(const DTenT &s, const size_t sdim);

template <typename TenElemT, typename QNT>
std::pair<size_t,size_t> CheckAndUpdateBoundaryMPSTensors(FiniteMPS<TenElemT, QNT> &,
                                                            const std::string&,
                                                            const size_t);

inline size_t CountLines(std::string filename);

template <typename TenElemT, typename QNT>
double TwoSiteFiniteVMPSSweep2(//also a overload
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const TwoSiteVMPSSweepParams &sweep_params,
    const size_t left_boundary,
    const size_t right_boundary,
    double& noise_start
);

template <typename QNT>
bool IsQNCovered(const QNSectorVec<QNT>&, const QNSectorVec<QNT>&);

/**
 Function to perform two-site noised update finite vMPS algorithm.

 @note The input MPS will be considered an empty one.
 @note The canonical center of MPS should be set at around left side
*/
template <typename TenElemT, typename QNT>
GQTEN_Double TwoSiteFiniteVMPS2( //same function name, overload by class of SweepParams 
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    TwoSiteVMPSSweepParams &sweep_params
){
    assert(mps.size() == mpo.size());
    std::cout << "***** Two-Site Noised Update VMPS Program (Single Processor) *****" << "\n";
    auto [left_boundary, right_boundary]=FiniteVMPSInit(mps,mpo,(SweepParams)sweep_params);
    std::cout << "Preseted noises: \t[";
    for(size_t i = 0; i < sweep_params.noises.size(); i++){
      std::cout << sweep_params.noises[i];
      if (i!=sweep_params.noises.size()-1) {
        std::cout << ", ";
      } else {
        std::cout << "]" << std::endl;
      }
    }
    GQTEN_Double e0;

    std::string file = "nf.json";
    std::ofstream ofs(file);

    ofs << "[\n";
    ofs.close();

    if (sweep_params.noises.size() == 0) { sweep_params.noises.push_back(0.0); }
    double noise_start;
    mps.LoadTen(left_boundary, GenMPSTenName(sweep_params.mps_path, left_boundary));
    mps.LoadTen(left_boundary+1, GenMPSTenName(sweep_params.mps_path, left_boundary+1));
    for (size_t sweep = 1; sweep <= sweep_params.sweeps; ++sweep) {
      if ((sweep - 1) < sweep_params.noises.size()) {
        noise_start = sweep_params.noises[sweep-1];
      }
      std::cout << "sweep " << sweep << std::endl;
      Timer sweep_timer("sweep");
      e0 = TwoSiteFiniteVMPSSweep2(mps, mpo, sweep_params, 
                                left_boundary, right_boundary, noise_start);
      sweep_timer.PrintElapsed();
      std::cout << std::endl;
    }
    mps.DumpTen(left_boundary, GenMPSTenName(sweep_params.mps_path, left_boundary), true);
    mps.DumpTen(left_boundary+1, GenMPSTenName(sweep_params.mps_path, left_boundary+1), true);
    
    ofs.open(file, std::ios_base::app);
    ofs << "]";
    ofs.close();
    return e0;
}


/**
Two-site (noised) update DMRG algorithm refer to 10.1103/PhysRevB.91.155115
*/
template <typename TenElemT, typename QNT>
double TwoSiteFiniteVMPSSweep2(//also a overload
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const TwoSiteVMPSSweepParams &sweep_params,
    const size_t left_boundary,
    const size_t right_boundary,
    double& noise_start
) {
  auto N = mps.size();
  using TenT = GQTensor<TenElemT, QNT>;
  TenVec<TenT> lenvs(N), renvs(N);
  double e0;

  double& noise_running = noise_start;
  for (size_t i = left_boundary; i < right_boundary-1; ++i) {
    //The last two site [right_boudary-1, right_boundary] will update when sweep back
    LoadRelatedTensTwoSiteAlg(mps, lenvs, renvs, i, 'r', sweep_params, left_boundary);    // note: here we need mps[i](do not need load),
                                                                            // mps[i+1](do not need load), mps[i+2](need load)
                                                                            // lenvs[i](do not need load), and mps[i+1]'s renvs
                                                                            // mps[i+1]'s renvs file can be removed
    e0 = TwoSiteFiniteVMPSUpdate2(
             mps,
             lenvs, renvs,
             mpo,
             sweep_params, 'r', i,
             noise_running
         );
    DumpRelatedTensTwoSiteAlg(mps, lenvs, renvs, i, 'r', sweep_params);    // note: here we need dump mps[i](free memory),
                                                                              // lenvs[i+1](without free memory)
  }

  for (size_t i = right_boundary; i > left_boundary+1; --i) {
    LoadRelatedTensTwoSiteAlg(mps, lenvs, renvs, i, 'l', sweep_params, right_boundary);
    e0 = TwoSiteFiniteVMPSUpdate2(
             mps,
             lenvs, renvs,
             mpo,
             sweep_params, 'l', i,
             noise_running
         );
    DumpRelatedTensTwoSiteAlg(mps, lenvs, renvs, i, 'l', sweep_params);
  }
  return e0;
}


/**  Single step for two site noised update.
This function includes below procedure:
- update `mps[target]` and `mps[next_site]` tensors according corresponding environment tensors and the mpo tensor, using lanczos algorithm;
- expand `mps[target]*mps[next_site]` and `mps[next_next_site]` by noise, if need;
- canonicalize mps to `mps[next_site]` by SVD, while truncate tensor `mps[target]` if need;
- generate the next environment in the direction.

When using this function, one must make sure memory at least contains `mps[target]` tensor,
`mps[next_site]`, its environment tensors, and `mps[next_next_site]`
*/
template <typename TenElemT, typename QNT>
double TwoSiteFiniteVMPSUpdate2(
    FiniteMPS<TenElemT, QNT> &mps,
    TenVec<GQTensor<TenElemT, QNT>> &lenvs,
    TenVec<GQTensor<TenElemT, QNT>> &renvs,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const TwoSiteVMPSSweepParams &sweep_params,
    const char dir,
    const size_t target_site,
    const double preset_noise
) {
  Timer update_timer("two_site_fvmps_update");

#ifdef GQMPS2_TIMING_MODE
  Timer preprocessing_timer("two_site_fvmps_initial_state");
#endif
  double noise = preset_noise;
  // Assign some parameters
  auto N = mps.size();
  std::vector<std::vector<size_t>> init_state_ctrct_axes;
  size_t svd_ldims;
  size_t lsite_idx, rsite_idx;
  size_t lenv_len, renv_len;
  std::string lblock_file, rblock_file;
  init_state_ctrct_axes = {{2}, {0}};
  svd_ldims = 2;
  switch (dir) {
    case 'r':
      lsite_idx = target_site;
      rsite_idx = target_site + 1;
      lenv_len = target_site;
      renv_len = N - (target_site + 2);
      break;
    case 'l':
      lsite_idx = target_site - 1;
      rsite_idx = target_site;
      lenv_len = target_site - 1;
      renv_len = N - target_site - 1;
      break;
    default:
      std::cout << "dir must be 'r' or 'l', but " << dir << std::endl;
      exit(1);
  }
 

  using TenT = GQTensor<TenElemT, QNT>;

  std::vector<TenT *>eff_ham(4);
  eff_ham[0] = lenvs(lenv_len);
  // Safe const casts for MPO local tensors.
  eff_ham[1] = const_cast<TenT *>(&mpo[lsite_idx]);
  eff_ham[2] = const_cast<TenT *>(&mpo[rsite_idx]);
  eff_ham[3] = renvs(renv_len);
  auto init_state = new TenT;
  Contract(&mps[lsite_idx], &mps[rsite_idx], init_state_ctrct_axes, init_state);
#ifdef GQMPS2_TIMING_MODE
   preprocessing_timer.PrintElapsed();
#endif
    // Lanczos
  Timer lancz_timer("two_site_fvmps_lancz");
  auto lancz_res = LanczosSolver(
                       eff_ham, init_state,
                       &eff_ham_mul_two_site_state,
                       sweep_params.lancz_params
                   );//Note here init_state is deleted
#ifdef GQMPS2_TIMING_MODE
  auto lancz_elapsed_time = lancz_timer.PrintElapsed();
#else
  auto lancz_elapsed_time = lancz_timer.Elapsed();
#endif


#ifdef GQMPS2_TIMING_MODE
  Timer expand_timer("two_site_fvmps_add_noise");
#endif

  bool need_expand(true);
  if (::fabs(noise) < 1e-10) {
    noise = 0.0;    //just for output
    need_expand = false;
  }else{
    const size_t physical_dim_l = mps[lsite_idx].GetShape()[1];
    const size_t physical_dim_r = mps[rsite_idx].GetShape()[1];
    const QNSectorVec<QNT>* qnscts_right;
    const QNSectorVec<QNT>* qnscts_left;
    Index<QNT> fused_index1, fused_index2;
    if (physical_dim_l == 2){
      qnscts_left  = &(mps[lsite_idx].GetIndexes()[0].GetQNScts());
    }else{
      std::vector<gqten::QNSctsOffsetInfo> qnscts_offset_info_list;
      fused_index1 = FuseTwoIndexAndRecordInfo(
            mps[lsite_idx].GetIndexes()[0],
            InverseIndex(mps[lsite_idx].GetIndexes()[1]),
            qnscts_offset_info_list
            );
      qnscts_left = &(fused_index1.GetQNScts());
    }

    if (physical_dim_r == 2){
      qnscts_right = &(mps[rsite_idx].GetIndexes()[2].GetQNScts());
    }else{
      std::vector<gqten::QNSctsOffsetInfo> qnscts_offset_info_list;
      fused_index2 = FuseTwoIndexAndRecordInfo(
            mps[rsite_idx].GetIndexes()[1],
            mps[rsite_idx].GetIndexes()[2],
            qnscts_offset_info_list
            );
      qnscts_right = &(fused_index2.GetQNScts());
    }

    if( dir == 'r' && 
        IsQNCovered(*qnscts_right, *qnscts_left) 
      ){
      noise = 0.0;            
      need_expand= false;
    }else if(dir == 'l' && 
        IsQNCovered(*qnscts_left, *qnscts_right)
      ){
      noise = 0.0;            
      need_expand= false;
    }
  }
  
  if (need_expand) {
    TwoSiteFiniteVMPSExpand(
        mps,
        lancz_res.gs_vec,
        eff_ham,
        dir,
        target_site,
        noise
    );
  }

#ifdef GQMPS2_TIMING_MODE
  expand_timer.PrintElapsed();
#endif


  // SVD and measure entanglement entropy
#ifdef GQMPS2_TIMING_MODE
  Timer svd_timer("two_site_fvmps_svd");
#endif

  TenT u, vt;
  using DTenT = GQTensor<GQTEN_Double, QNT>;
  DTenT s;
  GQTEN_Double actual_trunc_err;
  size_t D;
  SVD(
      lancz_res.gs_vec,
      svd_ldims, Div(mps[lsite_idx]),
      sweep_params.trunc_err, sweep_params.Dmin, sweep_params.Dmax,
      &u, &s, &vt, &actual_trunc_err, &D
  );
  delete lancz_res.gs_vec;
  auto ee = MeasureEE(s, D);

#ifdef GQMPS2_TIMING_MODE
  svd_timer.PrintElapsed();
#endif


  // Update MPS local tensor
#ifdef GQMPS2_TIMING_MODE
  Timer update_mps_ten_timer("two_site_fvmps_update_mps_ten");
#endif

  TenT the_other_mps_ten;
  switch (dir) {
    case 'r':
      mps[lsite_idx] = std::move(u);
      Contract(&s, &vt, {{1}, {0}}, &the_other_mps_ten);
      mps[rsite_idx] = std::move(the_other_mps_ten);
      break;
    case 'l':
      Contract(&u, &s, {{2}, {0}}, &the_other_mps_ten);
      mps[lsite_idx] = std::move(the_other_mps_ten);
      mps[rsite_idx] = std::move(vt);
      break;
    default:
      assert(false);
  }

#ifdef GQMPS2_TIMING_MODE
  update_mps_ten_timer.PrintElapsed();
#endif

// Update environment tensors
#ifdef GQMPS2_TIMING_MODE
  Timer update_env_ten_timer("two_site_fvmps_update_env_ten");
#endif

  switch (dir) {
    case 'r':{
      lenvs[lenv_len + 1] = std::move(UpdateSiteLenvs(lenvs[lenv_len], mps[target_site], mpo[target_site]));
    }break;
    case 'l':{
      renvs[renv_len + 1] = std::move(UpdateSiteRenvs(renvs[renv_len], mps[target_site], mpo[target_site])  );
    }break;
    default:
      assert(false);
  }


#ifdef GQMPS2_TIMING_MODE
  update_env_ten_timer.PrintElapsed();
#endif

  static TenT nf = TenT();
  static bool is_nf_initialized = false;
  if(!is_nf_initialized){
    mps.LoadTen(0, GenMPSTenName(sweep_params.mps_path, 0));
    Index<QNT> index_out_fermion = mps[0].GetIndexes()[1];
    Index<QNT> index_in_fermion = InverseIndex(index_out_fermion);
    nf = Tensor({index_in_fermion, index_out_fermion});
    nf({0,0}) = 2;
    nf({1,1}) = 1;
    nf({2,2}) = 1;
    nf({3,3}) = 0;
    mps.dealloc(0);
    is_nf_initialized = true;
  }
  const std::string file = "nf.json";

  size_t next_site = lsite_idx+rsite_idx-target_site;
  auto mps_ten_shape = mps[next_site].GetShape();
  if(mps_ten_shape[1]==4){
  //measure and dump singluar particle number
#ifdef GQMPS2_TIMING_MODE
  Timer measure_timer("single_site_fvmps_measure");
#endif
  MeasuResElem<TenElemT> nf_res = OneSiteOpAvg(
                        mps[next_site], nf,next_site,N );
  std::ofstream ofs(file, std::ios_base::app);
  ofs << "  [";
  DumpSites(ofs, nf_res.sites); DumpAvgVal(ofs, nf_res.avg);
  ofs << "],\n";
  ofs.close();
#ifdef GQMPS2_TIMING_MODE
  measure_timer.PrintElapsed();
#endif
  }

#ifdef GQMPS2_TIMING_MODE
  Timer max_sv_timer("single_site_fvmps_find_and_dump_max_sv");
#endif
  unsigned bond = (dir=='r')? (target_site+1):target_site;
  std::stringstream sstream;
  sstream << bond;
  std::string file_basename = sweep_params.mps_path+ "/sv_bond" + sstream.str();

  int pre_qnsct_num = CountLines(file_basename+".json") - 2;
  
  MaxSVInEachBlock( s, file_basename);
#ifdef GQMPS2_TIMING_MODE
  max_sv_timer.PrintElapsed();
#endif
  Index<QNT> middle_index = mps[lsite_idx].GetIndexes()[2];
  int post_qnsct_num = middle_index.GetQNSctNum();

  auto update_elapsed_time = update_timer.Elapsed();
  std::cout << "Site " << std::setw(4) << target_site
            << " E0 = " << std::setw(20) << std::setprecision(kLanczEnergyOutputPrecision) << std::fixed << lancz_res.gs_eng
            << " noise = " <<  std::setprecision(2) << std::scientific  << noise << std::fixed
            << " TruncErr = " << std::setprecision(2) << std::scientific << actual_trunc_err << std::fixed
            << " D = " << std::setw(5) << D
            << " Iter = " << std::setw(3) << lancz_res.iters
            << " LanczT = " << std::setw(8) << lancz_elapsed_time
            << " TotT = " << std::setw(8) << update_elapsed_time
            << " QnSctNumChange" << std::setw(3) << post_qnsct_num - pre_qnsct_num
            << " S = " << std::setw(10) << std::setprecision(7) << ee;
  std::cout << std::scientific << std::endl;
  return lancz_res.gs_eng;
}

/**
 * Counting how many lines in a file
 * if no file, return 0
 */
inline size_t CountLines(std::string filename){
  std::ifstream ReadFile(filename,std::ios::in);//read only
  if(ReadFile.fail()){
    return 0;
  }
  size_t n=0;
  std::string temp;

  while(getline(ReadFile,temp))
  {
    n++;
  }
  ReadFile.close();
  return n;
}

}//gqmps2
