#pragma once

#include "gqmps2/one_dim_tn/mps/finite_mps/finite_mps.h"    // FiniteMPS
#include "gqmps2/one_dim_tn/mps/finite_mps/finite_mps_measu.h"
#include "gqten/gqten.h"

#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace gqmps2 {
using namespace gqten;


// Measure one-site operator with specifying sites
/**
Measure a single one-site operator on specific sites of the finite MPS.

@tparam TenElemT Type of the tensor element.
@tparam QNT Quantum number type.

@param mps To-be-measured MPS.
@param op The single one-site operator.
@param sites The sites will be measured.
@param res_file_basename The basename of the output file.
*/
template <typename TenElemT, typename QNT>
MeasuRes<TenElemT> MeasureOneSiteOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const GQTensor<TenElemT, QNT> &op,
    const std::vector<size_t>& sites,
    const std::string &res_file_basename
) {
  size_t N = mps.size();
  size_t res_num = sites.size();
  MeasuRes<TenElemT> measu_res;
  measu_res.reserve(res_num);
  for (auto i: sites) {
    mps.Centralize(i);
    measu_res.push_back( OneSiteOpAvg(mps[i], op, i, N) );
  }
  DumpMeasuRes(measu_res, res_file_basename);
  return measu_res;
}



/**
Measure a list of one-site operators on specified sites of the finite MPS.

@tparam TenElemT Type of the tensor element.
@tparam QNT Quantum number type.

@param mps To-be-measured MPS.
@param ops A list of one-site operators.
@param sites The sites will be measured.
@param res_file_basename The basename of the output file.
*/
template <typename TenElemT, typename QNT>
MeasuResSet<TenElemT> MeasureOneSiteOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::vector<GQTensor<TenElemT, QNT>> &ops,
    const std::vector<size_t>& sites,
    const std::vector<std::string> &res_file_basenames
) {
  auto op_num = ops.size();
  assert(op_num == res_file_basenames.size());
  auto N = mps.size();
  size_t res_num = sites.size();
  MeasuResSet<TenElemT> measu_res_set(op_num);
  for (MeasuRes<TenElemT> &measu_res : measu_res_set) {
    measu_res.reserve(res_num);
  }
  for (size_t i: sites) {
    mps.Centralize(i);
    for (size_t j = 0; j < op_num; ++j) {
      measu_res_set[j].push_back( OneSiteOpAvg(mps[i], ops[j], i, N) ); 
    }
  }
  for (size_t i = 0; i < op_num; ++i) {
    DumpMeasuRes(measu_res_set[i], res_file_basenames[i]);
  }
  return measu_res_set;
}






/**
Measure a two-site operator without insertion operator.

@tparam TenElemT Type of the tensor element, real or complex.
@param mps To-be-measured MPS.
@param sites_set The indexes of the two physical operators with ascending order
       for each measure event. Its size defines the number of measure events.
@param res_file_basename The basename of the output file.
*/
template <typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureTwoSiteOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const GQTensor<TenElemT, QNT>& phys_ops1,
    const GQTensor<TenElemT, QNT>& phys_ops2,
    const std::vector<std::vector<size_t>> &sites_set,
    const std::string &res_file_basename
) {
  size_t measu_event_num = sites_set.size();
  for(auto& sites: sites_set){
    assert(sites.size()==2);
  };

  MeasuRes<TenElemT> measu_res(measu_event_num);
  for (std::size_t i = 0; i < measu_event_num; ++i) {
    size_t idx1 = sites_set[i][0];
    size_t idx2 = sites_set[i][1];
    mps.Centralize(idx1);
    measu_res[i] = TwoSiteOpAvg(mps, phys_ops1, phys_ops2, idx1, idx2);
  }
  DumpMeasuRes(measu_res, res_file_basename);
  return measu_res;
}

//two site operator average, No insertion case
template <typename TenElemT, typename QNT>
inline MeasuResElem<TenElemT> TwoSiteOpAvg(
    FiniteMPS<TenElemT, QNT> &mps,
    const GQTensor<TenElemT, QNT>& phys_op1,
    const GQTensor<TenElemT, QNT>& phys_op2,
    const size_t idx1,
    const size_t idx2
) {
  auto id_op_set = mps.GetSitesInfo().id_ops;
  std::vector<GQTensor<TenElemT, QNT>> ops(idx2-idx1+1);
  ops.front()=phys_op1;
  ops.back()=phys_op2;
  for(size_t i = 1;i<ops.size()-1;i++){
    ops[i] = id_op_set[idx1+i];
  }
  assert(idx1 < idx2);
  auto head_site = idx1;
  auto tail_site = head_site + ops.size() - 1;
  assert(tail_site == idx2);
  auto avg = OpsVecAvg(mps, ops, head_site, tail_site);

  return MeasuResElem<TenElemT>({idx1,idx2}, avg);
}

/**
Measure 4 point function, with some restriction.
This is a specially designed function used to calculate 4 point fermion correlation function 
in electron-phonon system, especially for d-wave/t-wave pair correlation.
@note The input phys_ops must be fermion operators (with matrix form of corresponding hardcore boson).
@note The sites_set can be devided into `Ly` groups. We have assumed `Ly=4` here. First two sites in every group is same.

@param mps To-be-measured MPS.
@param sites_set The indexes of the four physical operators with ascending order
       for each measure event. Its size defines the number of measure events.
@param res_file_basename The basename of the output file.

*/
template <typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureElectronPhonon4PointFunction(
  FiniteMPS<TenElemT, QNT> &mps,
  const std::vector<GQTensor<TenElemT, QNT>>& phys_ops,
  const std::vector<std::vector<size_t>> &sites_set,
  const size_t pesudosite_num,
  const std::string &res_file_basename
){
  const size_t Ly = 4; //change here when Ly of system changes
  std::cout <<  "note: Ly = " << Ly <<std::endl;
  assert(phys_ops.size() == 4);
  assert(sites_set.size()%Ly == 0);
  for(size_t i = 0; i < sites_set.size();i++){
    assert(sites_set[i].size() == 4);
  }
  const size_t measure_event_num = sites_set.size();
  const size_t measure_event_per_group_num = sites_set.size()/Ly;
  for(size_t group=0;group<Ly;group++){
    std::vector<std::vector<size_t>> sites_set_group(
      sites_set.cbegin() + group * measure_event_per_group_num,
      sites_set.cbegin() + (group+1) * measure_event_per_group_num 
      );
    size_t site1 = sites_set_group[0][0];
    size_t site2 = sites_set_group[0][1];
    assert(sites_set_group.size() == measure_event_per_group_num);
    for(auto iter = sites_set_group.begin(); iter< sites_set_group.cend();iter++){
      assert( (*iter)[0] = site1 );
      assert( (*iter)[1] = site2 );
      //TODO: check order: a. site1<site2<3<4 b. for every 3,4, ascending order
    }
  }
  MeasuRes<TenElemT> measure_res;
  measure_res.reserve(measure_event_num);
  for(size_t group=0; group<Ly;group++){
    std::vector<std::vector<size_t>> sites_set_group(
      sites_set.cbegin() + group * measure_event_per_group_num,
      sites_set.cbegin() + (group+1) * measure_event_per_group_num 
      );
      auto measure_res_group = MeasureElectronPhonon4PointFunctionGroup(mps,
        phys_ops,
        sites_set_group,
        pesudosite_num
        );
      measure_res.insert(measure_res.end(), measure_res_group.begin(),measure_res_group.end());
  }
  DumpMeasuRes(measure_res, res_file_basename);
  return measure_res;
}


template <typename TenElemT, typename QNT>
struct TempTensorWithContract2Ops{
  using TenT = GQTensor<TenElemT, QNT>;
  size_t idx; // The last site having be contract
  TenT* tmp_tensor;

  TempTensorWithContract2Ops(const size_t& idx, TenT* const & tmp_tensor):
  idx(idx), tmp_tensor(tmp_tensor) {}

  void MoveOnTo(
    const FiniteMPS<TenElemT, QNT> &mps,
    size_t site_to
  ){
    assert(site_to>=idx);
    for(idx=idx+1;idx<=site_to;idx++){
      CtrctMidTen(mps, idx, TenT(),TenT(),tmp_tensor);
    }
    idx--;
  }
};

template <typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureElectronPhonon4PointFunctionGroup(
  FiniteMPS<TenElemT, QNT> &mps,
  const std::vector<GQTensor<TenElemT, QNT>>& phys_ops,
  const std::vector<std::vector<size_t>> &sites_set,
  const size_t pesudosite_num//usually called Np in hxwang's program
){
  const size_t site1=sites_set[0][0];
  const size_t site2=sites_set[0][1];
  const size_t fermion_op_dim = mps[site1].GetIndexes()[1].dim();// should be 4 for SSH-Hubbard
  assert(pesudosite_num>0);
  const size_t bonson_op_dim = mps[site1+1].GetIndexes()[1].dim();// usually should be 2;
  assert(fermion_op_dim != bonson_op_dim);//we use dimension to differentiate boson site or fermion site 
  using Tensor = GQTensor<TenElemT, QNT>;
  using IndexT = Index<QNT>;
  Index<QNT> index_out_fermion = mps[0].GetIndexes()[1];
  Index<QNT> index_in_fermion = InverseIndex(index_out_fermion);
  static Tensor f = Tensor({index_in_fermion, index_out_fermion});
  static bool is_f_initial = false;
  if(!is_f_initial){
    f({0,0}) = 1;
    f({1,1}) = -1;
    f({2,2}) = -1;
    f({3,3}) = 1;
    is_f_initial = true;
  }

  mps.Centralize(site1);
  auto id_op_set = mps.GetSitesInfo().id_ops;
  //Contract on site1
  std::vector<size_t> head_mps_ten_ctrct_axes1 {1};
  std::vector<size_t> head_mps_ten_ctrct_axes2 {0, 2};
  std::vector<size_t> head_mps_ten_ctrct_axes3 {0, 1};
  GQTensor<TenElemT, QNT> temp_ten0;
  auto ptemp_ten = new GQTensor<TenElemT, QNT>;//delete when first called MoveOnTo, around line 283.
  Contract(
      &mps[site1], &phys_ops[0],
      {{1}, {0}},
      &temp_ten0
  );
  auto mps_ten_dag = Dag(mps[site1]);
  Contract(
      &temp_ten0, &mps_ten_dag,
      {head_mps_ten_ctrct_axes2, head_mps_ten_ctrct_axes3},
      ptemp_ten
  );

  for (size_t i = site1 + 1; i < site2; ++i) {
    if(mps[i].GetIndexes()[1].dim()==bonson_op_dim ){
      CtrctMidTen(mps, i, id_op_set[i], id_op_set[i], ptemp_ten);
    }else{
      CtrctMidTen(mps, i, f, id_op_set[i], ptemp_ten);
    }
  }
  CtrctMidTen(mps, site2, phys_ops[1], id_op_set[site2],ptemp_ten);
  TempTensorWithContract2Ops<TenElemT,QNT> tmp_tensor_with_c2o(site2, ptemp_ten);


  MeasuRes<TenElemT> measure_res(sites_set.size());
  for(size_t event = 0; event < sites_set.size() ;event++){//event means measure event
    size_t site3 = sites_set[event][2];
    size_t site4 = sites_set[event][3];
    tmp_tensor_with_c2o.MoveOnTo(mps,site3-1);//When first move on, ptemp_ten was deleted.
    ptemp_ten = new Tensor(*tmp_tensor_with_c2o.tmp_tensor);// deep copy
    CtrctMidTen(mps, site3, phys_ops[2], id_op_set[site3],ptemp_ten);
    for (size_t i = site3 + 1; i < site4; ++i) {
      if(mps[i].GetIndexes()[1].dim()==bonson_op_dim ){
        CtrctMidTen(mps, i, id_op_set[i], id_op_set[i], ptemp_ten);
      }else{
        CtrctMidTen(mps, i, f, id_op_set[i], ptemp_ten);
      }
    }
    // Deal with tail tensor.
    std::vector<size_t> tail_mps_ten_ctrct_axes1 {0, 1, 2};
    std::vector<size_t> tail_mps_ten_ctrct_axes2 {2, 0, 1};
    GQTensor<TenElemT, QNT> temp_ten2, temp_ten3, res_ten;
    Contract(&mps[site4], ptemp_ten, {{0}, {0}}, &temp_ten2);
    delete ptemp_ten;
    Contract(&temp_ten2, &phys_ops.back(), {{0}, {0}}, &temp_ten3);
    mps_ten_dag = std::move(Dag(mps[site4]));
    Contract(
      &temp_ten3, &mps_ten_dag,
      {tail_mps_ten_ctrct_axes1, tail_mps_ten_ctrct_axes2},
      &res_ten
    );
    measure_res[event]=MeasuResElem<TenElemT>(sites_set[event], res_ten());
  }
  delete tmp_tensor_with_c2o.tmp_tensor;
  return measure_res;
}



}//gqmps2