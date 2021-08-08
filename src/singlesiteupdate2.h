// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang
* Creation Date: 2021-6-12
*
* Description: GraceQ/MPS2 project. Implementation details for single-site algorithm.
* The difference between this version with the standarad version in GQMPS2 is 
* this version save the maximum singluar values in different quantum number sectors,
* and qunatum block infos. Also measure electron numbers when update fermion sites.
*/

/**
@file single_site_update2.h
@brief single-site finite variational MPS algorithm, magic change
*/

#ifndef GQMPS2_ALGORITM_VMPS_ONE_SITE_UPDATE_FINITE_VMPS_IMPL2_H
#define GQMPS2_ALGORITM_VMPS_ONE_SITE_UPDATE_FINITE_VMPS_IMPL2_H

#include "gqmps2/algorithm/vmps/single_site_update_finite_vmps.h"   // SingleVMPSSweepParams
#include "gqmps2/algorithm/vmps/two_site_update_finite_vmps.h"      // helper functions
#include "gqmps2/one_dim_tn/mpo/mpo.h"                              // MPO
#include "gqmps2/one_dim_tn/mps/finite_mps/finite_mps.h"  
#include "gqmps2/one_dim_tn/mps/finite_mps/finite_mps_measu.h"          // FiniteMPS
#include "gqmps2/utilities.h"                                       // IsPathExist, CreatPath
#include "gqmps2/one_dim_tn/framework/ten_vec.h"                    // TenVec
#include "gqmps2/consts.h"
#include "gqten/gqten.h"
#include "gqten/utility/timer.h"                                    // Timer

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <stdio.h>    // remove
#ifdef Release
#define NDEBUG
#endif
#include <assert.h>

namespace gqmps2 {
  using namespace gqten;

  using U1U1QN = QN<U1QNVal, U1QNVal>;
/**
Max singular values in every block
*/
struct MaxSVElem {
  MaxSVElem(void) = default;
  MaxSVElem(const std::vector<long> &qn_vals, const GQTEN_Double max_sv, const long degeneracy) :
    qn_vals(qn_vals), max_sv(max_sv), degeneracy(degeneracy) {}

  std::vector<long> qn_vals;  ///< QN vals label block
  GQTEN_Double max_sv;         ///< average of the observation.
  long degeneracy;
};

using MaxSVSet = std::vector<MaxSVElem>;

inline  void DumpMaxSVSet(
    const MaxSVSet &res,
    const std::string &basename
);

inline void DumpSVBlk(std::ofstream &ofs, const std::vector<long> &sites);
inline MaxSVSet MaxSVInEachBlock(const GQTensor<GQTEN_Double, U1U1QN>& s,const std::string& file_basename);
/**
Function to perform single-site update finite vMPS algorithm.

@note The input MPS will be considered an empty one.
@note The canonical center of MPS should be set at site 0.
*/
template <typename TenElemT, typename QNT>
GQTEN_Double SingleSiteFiniteVMPS2(
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    SingleVMPSSweepParams &sweep_params
){
    assert(mps.size() == mpo.size());

    std::cout << std::endl;
    std::cout << "=====> Sweep Parameter <=====" << std::endl;
    std::cout << "MPS/MPO size: \t " << mpo.size() << std::endl;
    std::cout << "The number of sweep times: \t " << sweep_params.sweeps << std::endl;
    std::cout << "Bond dimension: \t " << sweep_params.Dmin << "/" << sweep_params.Dmax << std::endl;
    std::cout << "Cut off truncation error: \t " <<sweep_params.trunc_err << std::endl;
    std::cout << "Lanczos max iterations \t" <<sweep_params.lancz_params.max_iterations << std::endl;
    std::cout << "Preseted noises: \t[";
    for(size_t i = 0; i < sweep_params.noises.size(); i++){
      std::cout << sweep_params.noises[i];
      if (i!=sweep_params.noises.size()-1) {
        std::cout << ", ";
      } else {
        std::cout << "]" << std::endl;
      }
    }
    std::cout << "MPS path: \t" << sweep_params.mps_path << std::endl;
    std::cout << "Temp path: \t" << sweep_params.temp_path << std::endl;

    // If the runtime temporary directory does not exit, create it and initialize
    // the left/right environments
    if (!IsPathExist(sweep_params.temp_path)) {
      CreatPath(sweep_params.temp_path);
      InitEnvs(mps, mpo, sweep_params.mps_path, sweep_params.temp_path, 1);
      std::cout << "no exsiting path " <<sweep_params.temp_path
                << ", thus progress created it and generated environment tensors."
                << std::endl;
    } else {
      std::cout << "finded exsiting path "<<sweep_params.temp_path
                << ", thus progress will use the present environment tensors."
                << std::endl;
    }

    std::string file = "nf.json";
    std::ofstream ofs(file);

    ofs << "[\n";
    ofs.close();


    GQTEN_Double e0;

    if (sweep_params.noises.size() == 0) { sweep_params.noises.push_back(0.0); }
    double noise_start;
    mps.LoadTen(0, GenMPSTenName(sweep_params.mps_path, 0));
    for (size_t sweep = 1; sweep <= sweep_params.sweeps; ++sweep) {
      if ((sweep - 1) < sweep_params.noises.size()) {
        noise_start = sweep_params.noises[sweep-1];
      }
      std::cout << "sweep " << sweep << std::endl;
      Timer sweep_timer("sweep");
      e0 = SingleSiteFiniteVMPSSweep2(mps, mpo, sweep_params, noise_start);
      sweep_timer.PrintElapsed();
      std::cout << std::endl;
    }
    mps.DumpTen(0, GenMPSTenName(sweep_params.mps_path, 0), true);
    
    ofs.open(file, std::ios_base::app);
    ofs << "]";
    ofs.close();
    
    return e0;
}


/**
Single-site update DMRG algorithm refer to 10.1103/PhysRevB.91.155115
*/
template <typename TenElemT, typename QNT>
double SingleSiteFiniteVMPSSweep2(
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const SingleVMPSSweepParams &sweep_params,
    double& noise_start
) {
  auto N = mps.size();
  using TenT = GQTensor<TenElemT, QNT>;
  TenVec<TenT> lenvs(N), renvs(N);
  double e0(0.0), actual_e0(0.0), actual_laststep_e0(0.0);

  const double alpha = sweep_params.alpha;
  const double noise_decrease = sweep_params.noise_decrease;
  const double noise_increase = sweep_params.noise_increase;
  const double max_noise = sweep_params.max_noise;

  double& noise_running = noise_start;
  for (size_t i = 0; i < N - 1; ++i) {
    LoadRelatedTensSingleSiteAlg(mps, lenvs, renvs, i, 'r', sweep_params);    // note: here we need mps[i](do not need load),
                                                                              // mps[i+1], lenvs[i](do not need load), and mps[i]'s renvs
                                                                              // mps[i]'s renvs can be removed
    actual_e0 = CalEnergyEptSingleSite(mps, mpo,lenvs, renvs, i);
    if ((actual_e0 - e0) <= 0.0) {
      // expand and truncate let the energy lower or not change
      // this case is very rare, but include the boundary mps tensor case
      // so we do nothing now
    } else if ((actual_e0 - e0) >= alpha*fabs(actual_laststep_e0-e0)) {
      // below two case suppose actual_laststep_e0-laststep_e0>0, usually it is right
      noise_running = noise_running*noise_decrease;
    } else {
      noise_running = std::min(noise_running*noise_increase, max_noise);
    }
    e0 = SingleSiteFiniteVMPSUpdate2(
             mps,
             lenvs, renvs,
             mpo,
             sweep_params, 'r', i,
             noise_running
         );
    actual_laststep_e0 = actual_e0;
    DumpRelatedTensSingleSiteAlg(mps, lenvs, renvs, i, 'r', sweep_params);    // note: here we need dump mps[i](free memory),
                                                                              // lenvs[i+1](without free memory)
  }

  for (size_t i = N-1; i > 0; --i) {
    LoadRelatedTensSingleSiteAlg(mps, lenvs, renvs, i, 'l', sweep_params);
    actual_e0 = CalEnergyEptSingleSite(mps, mpo,lenvs, renvs, i);
    if ((actual_e0 - e0) <= 0.0) {
    } else if ((actual_e0 - e0) >= alpha*fabs(actual_laststep_e0 - e0)) {
      noise_running = noise_running * noise_decrease;
    } else {
      noise_running = std::min(noise_running * noise_increase, max_noise);
    }
    e0 = SingleSiteFiniteVMPSUpdate2(
             mps,
             lenvs, renvs,
             mpo,
             sweep_params, 'l', i,
             noise_running
         );
    actual_laststep_e0 = actual_e0;
    DumpRelatedTensSingleSiteAlg(mps, lenvs, renvs, i, 'l', sweep_params);
  }
  return e0;
}


/**  Single step for single site update.
This function includes below procedure:
- update `mps[target]` tensors according corresponding environment tensors and the mpo tensor, using lanczos algorithm;
- expand `mps[target]` and `mps[next_site]` by noise, if need;
- canonicalize mps to `mps[next_site]` by SVD, while truncate tensor `mps[target]` if need;
- generate the next environment in the direction.

When using this function, one must make sure memory at least contains `mps[target]` tensor,
its environment tensors and `mps[next_site]`.
*/
template <typename TenElemT, typename QNT>
double SingleSiteFiniteVMPSUpdate2(
    FiniteMPS<TenElemT, QNT> &mps,
    TenVec<GQTensor<TenElemT, QNT>> &lenvs,
    TenVec<GQTensor<TenElemT, QNT>> &renvs,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const SingleVMPSSweepParams &sweep_params,
    const char dir,
    const size_t target_site,
    const double preset_noise
) {
  Timer update_timer("single_site_fvmps_update");

  double noise = preset_noise;
  auto N = mps.size();
  size_t lenv_len = target_site;
  size_t renv_len = N - target_site - 1;
  size_t svd_ldims;
  size_t next_site;
  switch (dir) {
    case 'r':
      svd_ldims = 2;
      next_site = target_site + 1;
      break;
    case 'l':
      svd_ldims = 1;
      next_site = target_site - 1;
      break;
    default:
      std::cout << "dir must be 'r' or 'l', but " << dir << std::endl;
      exit(3);
  }

  using TenT = GQTensor<TenElemT, QNT>;

    //particle number operator
  static TenT nf = TenT();
  static bool is_nf_initialized = false;
  if(!is_nf_initialized){
    Index<QNT> index_out_fermion = mps[0].GetIndexes()[1];
    Index<QNT> index_in_fermion = InverseIndex(index_out_fermion);
    nf = Tensor({index_in_fermion, index_out_fermion});
    nf({0,0}) = 2;
    nf({1,1}) = 1;
    nf({2,2}) = 1;
    nf({3,3}) = 0;
    is_nf_initialized = true;
  }
  const std::string file = "nf.json";

  std::vector<TenT *> eff_ham(3);
  eff_ham[0] = lenvs(lenv_len);
  eff_ham[1] = const_cast<TenT *>(mpo(target_site));    // Safe const casts for MPO local tensors.
  eff_ham[2] = renvs(renv_len);

  auto mps_ten_shape = mps[target_site].GetShape();
  Timer lancz_timer("single_site_fvmps_lancz");
  LanczosRes<TenT> lancz_res = LanczosSolver(
                       eff_ham,
                       mps(target_site),
                       &eff_ham_mul_single_site_state,
                       sweep_params.lancz_params
                   );                                   //note here mps(target_site) are destroyed.
#ifdef GQMPS2_TIMING_MODE
  auto lancz_elapsed_time = lancz_timer.PrintElapsed();
#else
  auto lancz_elapsed_time = lancz_timer.Elapsed();
#endif


  if(mps_ten_shape[1]==4){
  //measure and dump singluar particle number
#ifdef GQMPS2_TIMING_MODE
  Timer measure_timer("single_site_fvmps_measure");
#endif
  MeasuResElem<TenElemT> nf_res = OneSiteOpAvg(*lancz_res.gs_vec, nf,target_site,N );
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
  Timer noise_timer("single_site_fvmps_add_noise");
#endif

  bool need_expand(true);
  if (fabs(noise) < 1e-15) {
    need_expand = false;
  } else if (
      (target_site < N/2 && mps_ten_shape[0]*mps_ten_shape[1]<=mps_ten_shape[2]) ||
      (target_site > N/2 && mps_ten_shape[2]*mps_ten_shape[1]<=mps_ten_shape[0])
  ) {
    noise = 0.0;            //just for output
    need_expand= false;
  }
  if (need_expand) {
    SingleSiteFiniteVMPSExpand(
        mps,
        lancz_res.gs_vec,
        eff_ham,
        dir,
        target_site,
        noise
    );
    delete lancz_res.gs_vec;
  } else {
    mps(target_site) = lancz_res.gs_vec;
  }

#ifdef GQMPS2_TIMING_MODE
  noise_timer.PrintElapsed();
#endif

#ifdef GQMPS2_TIMING_MODE
  Timer svd_timer("single_site_fvmps_svd");
#endif

  TenT u, vt;
  GQTensor<GQTEN_Double, QNT> s;
  GQTEN_Double actual_trunc_err;
  size_t D;
  auto zero_div = Div(mps[target_site]) - Div(mps[target_site]);
  auto div_left = (dir=='r' ? Div(mps[target_site]) : zero_div);
  SVD(
      mps(target_site),
      svd_ldims, div_left,
      sweep_params.trunc_err, sweep_params.Dmin, sweep_params.Dmax,
      &u, &s, &vt, &actual_trunc_err, &D
  );
  auto ee = MeasureEE(s, D);

#ifdef GQMPS2_TIMING_MODE
  svd_timer.PrintElapsed();
#endif
#ifdef GQMPS2_TIMING_MODE
  Timer max_sv_timer("single_site_fvmps_find_and_dump_max_sv");
#endif
  unsigned bond = (dir=='r')? (target_site+1):target_site;
  std::stringstream sstream;
  sstream << bond;
  std::string file_basename = sweep_params.mps_path+ "/sv_bond" + sstream.str();
  MaxSVInEachBlock( s, file_basename);
#ifdef GQMPS2_TIMING_MODE
  max_sv_timer.PrintElapsed();
#endif

#ifdef GQMPS2_TIMING_MODE
  Timer update_mps_ten_timer("single_site_fvmps_update_mps_ten");
#endif

  TenT* temp_ten1 = new TenT();
  TenT* temp_ten2 = new TenT();
  switch(dir){
    case 'r':
      mps[target_site] = std::move(u);
      Contract(&s, &vt, {{1}, {0}}, temp_ten1);
      Contract(temp_ten1, mps(next_site), {{1}, {0}}, temp_ten2);
      delete temp_ten1;
      delete mps(next_site);
      mps(next_site) = temp_ten2;
      break;
    case 'l':
      mps[target_site] = std::move(vt);
      Contract(&u, &s, {{1}, {0}}, temp_ten1);
      Contract(mps(next_site), temp_ten1, {{2}, {0}}, temp_ten2);
      delete temp_ten1;
      delete mps(next_site);
      mps(next_site) = temp_ten2;
      break;
  }

#ifdef GQMPS2_TIMING_MODE
  update_mps_ten_timer.PrintElapsed();
#endif

// Update environment tensors
#ifdef GQMPS2_TIMING_MODE
  Timer update_env_ten_timer("single_site_fvmps_update_env_ten");
#endif

  switch (dir) {
    case 'r':{
      TenT temp1, temp2, lenv_ten;
      Contract(&lenvs[lenv_len], &mps[target_site], {{0}, {0}}, &temp1);
      Contract(&temp1, &mpo[target_site], {{0, 2}, {0, 1}}, &temp2);
      auto mps_ten_dag = Dag(mps[target_site]);
      Contract(&temp2, &mps_ten_dag, {{0 ,2}, {0, 1}}, &lenv_ten);
      lenvs[lenv_len + 1] = std::move(lenv_ten);
    }break;
    case 'l':{
      TenT temp1, temp2, renv_ten;
      Contract(&mps[target_site], eff_ham[2], {{2}, {0}}, &temp1);
      Contract(&temp1, &mpo[target_site], {{1, 2}, {1, 3}}, &temp2);
      auto mps_ten_dag = Dag(mps[target_site]);
      Contract(&temp2, &mps_ten_dag, {{3, 1}, {1, 2}}, &renv_ten);
      renvs[renv_len + 1] = std::move(renv_ten);
    }break;
    default:
      assert(false);
  }

#ifdef GQMPS2_TIMING_MODE
  update_env_ten_timer.PrintElapsed();
#endif


  auto update_elapsed_time = update_timer.Elapsed();
  std::cout << "Site " << std::setw(4) << target_site
            << " E0 = " << std::setw(20) << std::setprecision(kLanczEnergyOutputPrecision) << std::fixed << lancz_res.gs_eng
            << " noise = " <<  std::setprecision(2) << std::scientific  << noise << std::fixed
            << " TruncErr = " << std::setprecision(2) << std::scientific << actual_trunc_err << std::fixed
            << " D = " << std::setw(5) << D
            << " Iter = " << std::setw(3) << lancz_res.iters
            << " LanczT = " << std::setw(8) << lancz_elapsed_time
            << " TotT = " << std::setw(8) << update_elapsed_time
            << " S = " << std::setw(10) << std::setprecision(7) << ee;
  std::cout << std::scientific << std::endl;
  return lancz_res.gs_eng;
}


inline MaxSVSet MaxSVInEachBlock(
  const GQTensor<GQTEN_Double, U1U1QN>& s,
  const std::string& file_basename
){
  assert(s.Rank()==2);
  assert(Div(s) ==  Div(s)-Div(s) );

  const auto& bsdt_s=s.GetBlkSparDataTen();
  const auto& bidbm_s = bsdt_s.GetBlkIdxDataBlkMap();
  const GQTEN_Double* raw_data_pr_s = bsdt_s.GetActualRawDataPtr();
  std::cout << "qn block num of singular value" <<bidbm_s.size()<<std::endl;
  MaxSVSet max_sv_set; 
  max_sv_set.reserve(bidbm_s.size());
  for(auto& [blkidx, datablk] : bidbm_s ){
    const U1U1QN qn = datablk.GetQNBlkInfo().qnscts[0].GetQn();
    const long degeneracy=datablk.shape[0];
    int Nval = qn.GetQNVal(0).GetVal();
    int Szval = qn.GetQNVal(1).GetVal();
    const GQTEN_Double* pmax_sv = raw_data_pr_s+datablk.data_offset;
    GQTEN_Double max_sv = *(pmax_sv);
    max_sv_set.push_back(MaxSVElem({Nval, Szval}, max_sv, degeneracy));
  }
  DumpMaxSVSet(max_sv_set, file_basename);
  return max_sv_set;
}



// Data dump.
inline void DumpMaxSVSet(
    const MaxSVSet &res,
    const std::string &basename
) {
  auto file = basename + ".json";
  std::ofstream ofs(file);

  ofs << "[\n";

  for (auto it = res.begin(); it != res.end(); ++it) {
    auto &measu_res_elem = *it;

    ofs << "  [";

    DumpSVBlk(ofs, measu_res_elem.qn_vals); 
    DumpAvgVal(ofs, measu_res_elem.max_sv);
    ofs <<",";
    DumpAvgVal(ofs, measu_res_elem.degeneracy);

    if (it == res.end()-1) {
      ofs << "]\n";
    } else {
      ofs << "],\n";
    }
  }

  ofs << "]";

  ofs.close();
}


inline void DumpSVBlk(std::ofstream &ofs, const std::vector<long> &sites) {
  ofs << "[";
  for (auto it = sites.begin(); it != sites.end()-1; ++it) {
    ofs << *it << ", ";
  }
  ofs << sites.back();
  ofs << "], ";
}

inline void DumpAvgVal(std::ofstream &ofs, const long num) {
  ofs << std::setw(14) << num;
}

} /* gqmps2 */
#endif //GQMPS2_ALGORITM_VMPS_ONE_SITE_UPDATE_FINITE_VMPS_IMPL_H
