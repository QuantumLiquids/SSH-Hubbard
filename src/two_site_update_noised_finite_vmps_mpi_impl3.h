// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2021-09-11
*
* Description: GraceQ/MPS2 project. Two-site update noised finite size vMPS with MPI Paralization
*
* continue version, only sweep = 1 time, from middle,
* The API different with standard MPS2 API by Add "2" at the end of functions.
*/

/**
@file two_site_update_noised_finite_vmps_mpi_impl.h
@brief Two-site update noised finite size vMPS with MPI Paralization
*/
#ifndef GQMPS2_ALGO_MPI_VMPS_TWO_SITE_UPDATE_NOISED_FINITE_VMPS_MPI_IMPL2_H
#define GQMPS2_ALGO_MPI_VMPS_TWO_SITE_UPDATE_NOISED_FINITE_VMPS_MPI_IMPL2_H



#include <cstdlib>
#include "gqten/gqten.h"
#include "gqmps2/algorithm/lanczos_solver.h"                        //LanczosParams
#include "boost/mpi.hpp"                                            //boost::mpi
#include "gqmps2/algo_mpi/framework.h"                              //VMPSORDER
#include "gqmps2/algo_mpi/vmps/vmps_mpi_init.h"                     //MPI vmps initial
#include "gqmps2/algo_mpi/vmps/two_site_update_finite_vmps_mpi.h"   //TwoSiteMPIVMPSSweepParams
#include "gqmps2/algo_mpi/vmps/two_site_update_noised_finite_vmps_mpi.h" //TwoSiteMPINoisedVMPSSweepParams
#include "gqmps2/algo_mpi/lanczos_solver_mpi.h"                     //MPI Lanczos solver
#include "gqmps2/algo_mpi/vmps/two_site_update_finite_vmps_mpi_impl.h" //SlaveTwoSiteFiniteVMPS
#include "gqmps2/algo_mpi/vmps/two_site_update_noised_finite_vmps_mpi_impl.h" //Load related tensors
#include <thread>                                                       //thread


namespace gqmps2 {
using namespace gqten;

//forward decelaration
template <typename TenElemT, typename QNT>
inline void LoadRelatedTensOnTwoSiteAlgWhenNoisedRightMoving(
    FiniteMPS<TenElemT, QNT> &mps,
    TenVec<GQTensor<TenElemT, QNT>> &lenvs,
    TenVec<GQTensor<TenElemT, QNT>> &renvs,
    const size_t target_site,
    const size_t left_boundary,
    const TwoSiteMPINoisedVMPSSweepParams &sweep_params
);


template <typename TenElemT, typename QNT>
inline void LoadRelatedTensOnTwoSiteAlgWhenNoisedLeftMoving(
    FiniteMPS<TenElemT, QNT> &mps,
    TenVec<GQTensor<TenElemT, QNT>> &lenvs,
    TenVec<GQTensor<TenElemT, QNT>> &renvs,
    const size_t target_site,
    const size_t right_boundary,
    const TwoSiteMPINoisedVMPSSweepParams &sweep_params
);


/**
 * @note It's better master the tensor manipulation thread number = slave's -2.
 *       For the other threads used to read/dump tensors.
 */
template <typename TenElemT, typename QNT>
inline GQTEN_Double TwoSiteFiniteVMPS2(
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    TwoSiteMPINoisedVMPSSweepParams &sweep_params,
    mpi::communicator& world,
    const size_t start_site,
    const char start_direction
){
  GQTEN_Double e0(0.0);
  if(world.rank()== kMasterRank){
    e0 = MasterTwoSiteFiniteVMPS2(mps,mpo,sweep_params,world, start_site, start_direction);
  }else{
    SlaveTwoSiteFiniteVMPS<TenElemT, QNT>(world);
  }
  return e0;
}


template <typename TenElemT, typename QNT>
GQTEN_Double MasterTwoSiteFiniteVMPS2(
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    TwoSiteMPINoisedVMPSSweepParams &sweep_params,
    mpi::communicator world,
    const size_t start_site,
    const char start_direction
) {
  sweep_params.sweeps = 1;
  std::cout << "Note program set sweep time = 1!" << std::endl;
  assert(mps.size() == mpo.size());
  std::cout << "***** Two-Site Noised Update VMPS FIX Program (with MPI Parallel) *****" << "\n";
  MasterBroadcastOrder(program_start, world );
  std::cout << "=====> Checking and updating boundary tensors =====>" << std::endl;
  auto [left_boundary, right_boundary] = CheckAndUpdateBoundaryMPSTensors(
      mps,
      sweep_params.mps_path,
      sweep_params.Dmax
  );
  using TenT = GQTensor<TenElemT, QNT>;
  UpdateBoundaryEnvs(mps, mpo, sweep_params.mps_path,
                     sweep_params.temp_path, left_boundary, right_boundary, 2 );
  ///< one more left environment is need to update
  std::string Tfile =GenEnvTenName("l", left_boundary, sweep_params.temp_path);
  TenT lenv;
  ReadGQTensorFromFile(lenv, Tfile);
  mps.LoadTen(left_boundary,
      GenMPSTenName(sweep_params.mps_path, left_boundary));
  lenv = std::move(UpdateSiteLenvs(lenv, mps[left_boundary], mpo[left_boundary]));
  Tfile =GenEnvTenName("l", left_boundary+1, sweep_params.temp_path);
  WriteGQTensorTOFile(lenv, Tfile);

  std::cout << "Preseted noises: \t[";
  for(size_t i = 0; i < std::min(sweep_params.sweeps, sweep_params.noises.size()); i++){
    std::cout << sweep_params.noises[i];
    if (i!=sweep_params.noises.size()-1) {
      std::cout << ", ";
    } else {
      std::cout << "]" << std::endl;
    }
  }

  std::string file = "nf.json";
  std::ofstream ofs(file);
  ofs << "[\n";
  ofs.close();
  if (sweep_params.noises.empty()) { sweep_params.noises.push_back(0.0); }
  double e0(0.0);
  double noise;
  
  for (size_t sweep = 1; sweep <= sweep_params.sweeps; ++sweep) {
    if ((sweep - 1) < sweep_params.noises.size()) {
      noise = sweep_params.noises[sweep-1];
    }
    std::cout << "sweep " << sweep << std::endl;
    Timer sweep_timer("sweep");
    if(start_direction == 'r' ){
      e0 = TwoSiteFiniteVMPSSweep2_StartToRight(mps, mpo, sweep_params,
                                                left_boundary, right_boundary,
                                                noise,  world, start_site);
    } else if(start_direction =='l'){
      e0 = TwoSiteFiniteVMPSSweep2_StartToLeft(mps, mpo, sweep_params,
                                                left_boundary, right_boundary,
                                                noise,  world, start_site);
    } else {
      std::cout << "start_direction = " << start_direction << std::endl;
      exit(1);
    }



    sweep_timer.PrintElapsed();
    std::cout << "\n";
  }
  mps.LeftCanonicalizeTen(left_boundary);
  mps.DumpTen(left_boundary, GenMPSTenName(sweep_params.mps_path, left_boundary), true);
  mps.DumpTen(left_boundary+1, GenMPSTenName(sweep_params.mps_path, left_boundary+1), true);
  ofs.open(file, std::ios::in);
  ofs.seekp(-2, std::ios::end);
  ofs << "\n]";
  ofs.close();
  MasterBroadcastOrder(program_final, world);
  return e0;
}

template <typename TenElemT, typename QNT>
double TwoSiteFiniteVMPSSweep2_StartToRight(
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const TwoSiteMPINoisedVMPSSweepParams &sweep_params,
    const size_t left_boundary,
    const size_t right_boundary,
    const double noise,
    mpi::communicator world,
    const size_t start_site
) {
  std::cout << "To right"<< std::endl;
  auto N = mps.size();
  TwoSiteMPIVMPSSweepParams sweep_params_no_noise = (TwoSiteMPIVMPSSweepParams) sweep_params;
  using TenT = GQTensor<TenElemT, QNT>;
  TenVec<TenT> lenvs(N - 1);
  TenVec<TenT> renvs(N - 1);
  double e0;
  const size_t update_site_size = right_boundary - left_boundary -1;
  std::thread load_related_tens_thread;
  std::thread dump_related_tens_thread;
  mps.LoadTen(start_site, GenMPSTenName(sweep_params.mps_path, start_site));
  mps.LoadTen(start_site+1, GenMPSTenName(sweep_params.mps_path, start_site+1));
  lenvs.LoadTen(start_site, GenEnvTenName("l", start_site, sweep_params.temp_path) );
  LoadRelatedTensOnTwoSiteAlgWhenNoisedRightMoving(mps, lenvs, renvs, start_site, left_boundary, sweep_params);
  for (size_t i = start_site; i <= right_boundary - 2; ++i) {
    // Load to-be-used tensors
    if( i <  right_boundary - 2 ){
      load_related_tens_thread = std::thread(
          LoadRelatedTensOnTwoSiteAlgWhenNoisedRightMoving<TenElemT, QNT>,
          std::ref(mps),
          std::ref(lenvs),
          std::ref(renvs),
          i + 1, //note here is different,
          left_boundary,
          std::ref(sweep_params)
      );
    }
    if(i==start_site){
      TwoSiteMPINoisedVMPSSweepParams sweep_params2 = sweep_params;
      sweep_params2.lancz_params.max_iterations = 100;
      e0 = MasterTwoSiteFiniteVMPSUpdate2(mps, lenvs, renvs, mpo, sweep_params2, 'r', i, noise,world);
    }else{
      e0 = MasterTwoSiteFiniteVMPSUpdate2(mps, lenvs, renvs, mpo, sweep_params, 'r', i, noise,world);
    }

    // Dump related tensor to HD and remove unused tensor from RAM
    if( i > start_site ){
      dump_related_tens_thread.join();
    }
    dump_related_tens_thread= std::thread(
        DumpRelatedTensOnTwoSiteAlgWhenRightMoving<TenElemT, QNT>,
        std::ref(mps),
        std::ref(lenvs),
        std::ref(renvs),
        i,
        std::ref(sweep_params_no_noise)
    );
    if ( i < right_boundary - 2 ){
      load_related_tens_thread.join();
    }
  }
  dump_related_tens_thread.join();

  LoadRelatedTensOnTwoSiteAlgWhenNoisedLeftMoving(mps, lenvs, renvs, right_boundary, right_boundary, sweep_params);
  for (size_t i = right_boundary; i >= left_boundary+2; --i) {
    if(i > left_boundary+2){
      load_related_tens_thread = std::thread(
          LoadRelatedTensOnTwoSiteAlgWhenNoisedLeftMoving<TenElemT, QNT>,
          std::ref(mps),
          std::ref(lenvs),
          std::ref(renvs),
          i - 1, //note here is different,
          right_boundary,
          std::ref(sweep_params)
      );
    }
    // LoadRelatedTensOnTwoSiteAlgWhenNoisedLeftMoving(mps, lenvs, renvs, i, right_boundary, sweep_params);
    e0 = MasterTwoSiteFiniteVMPSUpdate2(mps, lenvs, renvs, mpo, sweep_params, 'l', i, noise, world);
    if( i < right_boundary ) {
      dump_related_tens_thread.join();
    }
    dump_related_tens_thread = std::thread(
        DumpRelatedTensOnTwoSiteAlgWhenLeftMoving<TenElemT, QNT>,
        std::ref(mps),
        std::ref(lenvs),
        std::ref(renvs),
        i,
        std::ref(sweep_params_no_noise)
    );
    if(i > left_boundary+2){
      load_related_tens_thread.join();
    }
  }
  dump_related_tens_thread.join();
  return e0;
}




template <typename TenElemT, typename QNT>
double TwoSiteFiniteVMPSSweep2_StartToLeft(
    FiniteMPS<TenElemT, QNT> &mps,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const TwoSiteMPINoisedVMPSSweepParams &sweep_params,
    const size_t left_boundary,
    const size_t right_boundary,
    const double noise,
    mpi::communicator world,
    const size_t start_site
) {
  std::cout <<"To left."<<std::endl;
  auto N = mps.size();
  TwoSiteMPIVMPSSweepParams sweep_params_no_noise = (TwoSiteMPIVMPSSweepParams) sweep_params;
  using TenT = GQTensor<TenElemT, QNT>;
  TenVec<TenT> lenvs(N - 1);
  TenVec<TenT> renvs(N - 1);
  double e0;
  std::thread load_related_tens_thread;
  std::thread dump_related_tens_thread;
  mps.LoadTen(start_site - 1, GenMPSTenName(sweep_params.mps_path, start_site-1));
  mps.LoadTen(start_site, GenMPSTenName(sweep_params.mps_path, start_site));
  lenvs.LoadTen(start_site -1, GenEnvTenName("l", start_site -1, sweep_params.temp_path) );
  renvs.LoadTen(N - 1 - (start_site ), GenEnvTenName("r", (N-1)-(start_site), sweep_params.temp_path));

  LoadRelatedTensOnTwoSiteAlgWhenNoisedLeftMoving(mps, lenvs, renvs, start_site , right_boundary, sweep_params);
  std::cout << "mps[" << start_site-1 <<"] :" << std::endl;
  mps(start_site -1 )->ConciseShow();
//  mps(start_site -1)->GetIndexes()[0]->Show();


  std::cout << "mps[" << start_site <<"] :" << std::endl;
  mps(start_site )->ConciseShow();


  std::cout << "lenvs["<< start_site - 1 <<"] :" << std::endl;
  lenvs(start_site -1 )->ConciseShow();



  for (size_t i = start_site; i >= left_boundary+2; --i) {
    if(i > left_boundary+2){
      load_related_tens_thread = std::thread(
          LoadRelatedTensOnTwoSiteAlgWhenNoisedLeftMoving<TenElemT, QNT>,
          std::ref(mps),
          std::ref(lenvs),
          std::ref(renvs),
          i - 1, //note here is different,
          right_boundary,
          std::ref(sweep_params)
      );
    }
    if(i==start_site){
      TwoSiteMPINoisedVMPSSweepParams sweep_params2 = sweep_params;
      sweep_params2.lancz_params.max_iterations = 100;
      e0 = MasterTwoSiteFiniteVMPSUpdate2(mps, lenvs, renvs, mpo, sweep_params2, 'l', i, noise,world);
    }else{
      e0 = MasterTwoSiteFiniteVMPSUpdate2(mps, lenvs, renvs, mpo, sweep_params, 'l', i, noise,world);
    }
    if( i < start_site ) {
      dump_related_tens_thread.join();
    }
    dump_related_tens_thread = std::thread(
        DumpRelatedTensOnTwoSiteAlgWhenLeftMoving<TenElemT, QNT>,
        std::ref(mps),
        std::ref(lenvs),
        std::ref(renvs),
        i,
        std::ref(sweep_params_no_noise)
    );
    if(i > left_boundary+2){
      load_related_tens_thread.join();
    }
  }
  dump_related_tens_thread.join();
  return e0;
}
















template <typename TenElemT, typename QNT>
double MasterTwoSiteFiniteVMPSUpdate2(
    FiniteMPS<TenElemT, QNT> &mps,
    TenVec<GQTensor<TenElemT, QNT>> &lenvs,
    TenVec<GQTensor<TenElemT, QNT>> &renvs,
    const MPO<GQTensor<TenElemT, QNT>> &mpo,
    const TwoSiteMPINoisedVMPSSweepParams &sweep_params,
    const char dir,
    const size_t target_site,
    double noise,
    mpi::communicator& world
) {
  //master
  Timer update_timer("two_site_fvmps_update");
#ifdef GQMPS2_TIMING_MODE
  Timer initialize_timer("two_site_fvmps_setup_and_initial_state");
#endif
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

  // Lanczos
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
  initialize_timer.PrintElapsed();
#endif
  Timer lancz_timer("two_site_fvmps_lancz");
  MasterBroadcastOrder(lanczos, world);
  for(size_t i = 0; i < 4; i++) {
    std::cout << " raw data of eff_ham[" << i <<"] = " << eff_ham[i]->GetBlkSparDataTen().GetActualRawDataSize() << std::endl;
  }
  auto lancz_res = MasterLanczosSolver(
      eff_ham, init_state,
      sweep_params.lancz_params,
      world
  );
#ifdef GQMPS2_TIMING_MODE
  auto lancz_elapsed_time = lancz_timer.PrintElapsed();
#else
  auto lancz_elapsed_time = lancz_timer.Elapsed();
#endif

  bool need_expand(true);
  if (fabs(noise) < 1e-10) {
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
    if(dir=='r'){
      MasterBroadcastOrder(contract_for_right_moving_expansion, world);
      MasterTwoSiteFiniteVMPSRightMovingExpand(
          mps,
          lancz_res.gs_vec,
          eff_ham,
          target_site,
          noise,
          world
      );
    }else{
      MasterBroadcastOrder(contract_for_left_moving_expansion, world);
      MasterTwoSiteFiniteVMPSLeftMovingExpand(
          mps,
          lancz_res.gs_vec,
          eff_ham,
          target_site,
          noise,
          world
      );
    }
  }

  // SVD and measure entanglement entropy
#ifdef GQMPS2_TIMING_MODE
  Timer svd_timer("two_site_fvmps_svd");
#endif

  TenT u, vt;
  using DTenT = GQTensor<GQTEN_Double, QNT>;
  DTenT s;
  GQTEN_Double actual_trunc_err;
  size_t D;
  MasterBroadcastOrder(svd, world);
  MPISVDMaster(
      lancz_res.gs_vec,
      svd_ldims, Div(mps[lsite_idx]),
      sweep_params.trunc_err, sweep_params.Dmin, sweep_params.Dmax,
      &u, &s, &vt, &actual_trunc_err, &D,
      world
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
  switch (dir) {
    case 'r':{
      MasterBroadcastOrder(growing_left_env, world);
      lenvs(lenv_len + 1) = MasterGrowLeftEnvironment(lenvs[lenv_len], mpo[target_site],mps[target_site], world);
    }break;
    case 'l':{
      MasterBroadcastOrder(growing_right_env, world);
      renvs(renv_len + 1) = MasterGrowRightEnvironment(*eff_ham[3], mpo[target_site],mps[target_site], world);
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
            << " TruncErr = " << std::setprecision(2) << std::scientific << actual_trunc_err << std::fixed
            << " D = " << std::setw(5) << D
            << " Iter = " << std::setw(3) << lancz_res.iters
            << " LanczT = " << std::setw(8) << lancz_elapsed_time
            << " TotT = " << std::setw(8) << update_elapsed_time
            << " S = " << std::setw(10) << std::setprecision(7) << ee;
  std::cout << std::scientific << std::endl;
  return lancz_res.gs_eng;
}

}//gqmps2
#endif
