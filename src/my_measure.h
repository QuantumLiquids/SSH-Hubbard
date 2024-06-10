#ifndef MY_MEASURE_H
#define MY_MEASURE_H

#include "qlmps/one_dim_tn/mps/finite_mps/finite_mps.h"    // FiniteMPS
#include "qlmps/one_dim_tn/mps/finite_mps/finite_mps_measu.h"
#include "qlten/qlten.h"

#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "boost/mpi.hpp"

namespace qlmps {
using namespace qlten;

//helper
///< @note tensors from site 0 to left_boundary + 1 are loaded
template<typename TenElemT, typename QNT>
size_t FindLeftBoundary(FiniteMPS<TenElemT, QNT> &mps) {
  assert(mps.empty());
  size_t N = mps.size();
  const std::string mps_path = kMpsPath; //only for usual usage.
  const size_t left_middle_site = N / 2 - 1; //only for system large case, almost always work in phonon project
  using TenT = QLTensor<TenElemT, QNT>;
  mps.LoadTen(0, GenMPSTenName(mps_path, 0));
  size_t left_boundary(0);
  for (size_t i = 0; i < left_middle_site; i++) {
    mps.LoadTen(i + 1, GenMPSTenName(mps_path, i + 1));

    TenT &mps_ten = mps[i];
    ShapeT mps_ten_shape = mps_ten.GetShape();
    if (mps_ten_shape[0] * mps_ten_shape[1] > mps_ten_shape[2]) {
      left_boundary = i;
      size_t Dmax = mps_ten_shape[2];
      std::cout << "Bond dimension should be D = " << Dmax << "\n";
      std::cout << "left boundary site = " << left_boundary << "\n";
      break;
    }
    if (i == left_middle_site - 1) {
      left_boundary = i;
    }
  }
  return left_boundary;
}

/**
Measure a single one-site operator on specific sites of the finite MPS.
Memory are optimized and the input mps should be a empty mps.
The disk data will not change when and after measuring.

@tparam TenElemT Type of the tensor element.
@tparam QNT Quantum number type.

@param mps To-be-measured MPS.
@param op The single one-site operator.
@param sites The sites will be measured.
@param res_file_basename The basename of the output file.
*/
template<typename TenElemT, typename QNT>
MeasuRes<TenElemT> MeasureOneSiteOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const QLTensor<TenElemT, QNT> &op,
    const std::vector<size_t> &sites,
    const std::string &res_file_basename
) {
  size_t N = mps.size();
  size_t res_num = sites.size();
  MeasuRes<TenElemT> measu_res;
  measu_res.reserve(res_num);

  //Find the canonical center. We suppose the center = first site which is not complete orthogonal transformation + 1
  const size_t left_boundary = FindLeftBoundary(mps);
  const size_t initial_center = left_boundary + 1;
  //below we suppose sites[0] == 0
  for (size_t i = initial_center; i > 0; i--) {
    mps.RightCanonicalizeTen(i);
  }

  const std::string mps_path = kMpsPath;

  for (size_t i = 0; i < sites.size(); i++) {
    const size_t site = sites[i];
    if (i == 0) {
      measu_res.push_back(OneSiteOpAvg(mps[site], op, site, N));
      continue;
    }

    const size_t last_site = sites[i - 1];
    for (size_t j = last_site; j < site; j++) {
      if (j >= initial_center) {
        mps.LoadTen(j + 1, GenMPSTenName(mps_path, j + 1));
      }
      mps.LeftCanonicalizeTen(j);
      mps.dealloc(j);
    }
    measu_res.push_back(OneSiteOpAvg(mps[site], op, site, N));
    std::cout << "measured site " << site << "\n";
  }
  mps.dealloc(sites.back());
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
template<typename TenElemT, typename QNT>
MeasuResSet<TenElemT> MeasureOneSiteOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::vector<QLTensor<TenElemT, QNT>> &ops,
    const std::vector<size_t> &sites,
    const std::vector<std::string> &res_file_basenames
) {
  auto op_num = ops.size();
  assert(op_num == res_file_basenames.size());
  auto N = mps.size();
  size_t res_num = sites.size();
  MeasuResSet<TenElemT> measu_res_set(op_num);
  for (MeasuRes<TenElemT> &measu_res: measu_res_set) {
    measu_res.reserve(res_num);
  }

  //Find the canonical center. We suppose the center = first site which is not complete orthogonal transformation + 1
  const size_t left_boundary = FindLeftBoundary(mps);
  const size_t initial_center = left_boundary + 1;
  //below we suppose sites[0] == 0
  for (size_t i = initial_center; i > 0; i--) {
    mps.RightCanonicalizeTen(i);
  }
  const std::string mps_path = kMpsPath;

  for (size_t i = 0; i < sites.size(); i++) {
    const size_t site = sites[i];
    if (i == 0) {
      for (size_t j = 0; j < op_num; ++j) {
        measu_res_set[j].push_back(OneSiteOpAvg(mps[site], ops[j], site, N));
      }
      continue;
    }

    const size_t last_site = sites[i - 1];
    for (size_t j = last_site; j < site; j++) {
      if (j >= initial_center) {
        mps.LoadTen(j + 1, GenMPSTenName(mps_path, j + 1));
      }
      mps.LeftCanonicalizeTen(j);
      mps.dealloc(j);
    }
    for (size_t j = 0; j < op_num; ++j) {
      measu_res_set[j].push_back(OneSiteOpAvg(mps[site], ops[j], site, N));
    }
    std::cout << "measured site " << site << std::endl;
  }
  mps.dealloc(sites.back());
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
template<typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureTwoSiteOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const QLTensor<TenElemT, QNT> &phys_ops1,
    const QLTensor<TenElemT, QNT> &phys_ops2,
    const std::vector<std::vector<size_t>> &sites_set,
    const std::string &res_file_basename
) {
  size_t measu_event_num = sites_set.size();
  for (auto &sites: sites_set) {
    assert(sites.size() == 2);
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
template<typename TenElemT, typename QNT>
inline MeasuResElem<TenElemT> TwoSiteOpAvg(
    FiniteMPS<TenElemT, QNT> &mps,
    const QLTensor<TenElemT, QNT> &phys_op1,
    const QLTensor<TenElemT, QNT> &phys_op2,
    const size_t idx1,
    const size_t idx2
) {
  auto id_op_set = mps.GetSitesInfo().id_ops;
  std::vector<QLTensor<TenElemT, QNT>> ops(idx2 - idx1 + 1);
  ops.front() = phys_op1;
  ops.back() = phys_op2;
  for (size_t i = 1; i < ops.size() - 1; i++) {
    ops[i] = id_op_set[idx1 + i];
  }
  assert(idx1 < idx2);
  auto head_site = idx1;
  auto tail_site = head_site + ops.size() - 1;
  assert(tail_site == idx2);
  auto avg = OpsVecAvg(mps, ops, head_site, tail_site);

  return MeasuResElem<TenElemT>({idx1, idx2}, avg);
}

/**
Measure a two-site operator without insertion operator.
MPI version, memory optimized version

@tparam TenElemT Type of the tensor element, real or complex.
@param mps To-be-measured MPS.
@param sites_set The indexes of the two physical operators with ascending order
       for each measure event. Its size defines the number of measure events.
@param res_file_basename The basename of the output file.
*/
template<typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureTwoSiteOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const QLTensor<TenElemT, QNT> &phys_ops1,
    const QLTensor<TenElemT, QNT> &phys_ops2,
    const std::vector<std::vector<size_t>> &sites_set,
    const size_t Ly,
    const std::string &res_file_basename,
    const boost::mpi::communicator &world
) {
  assert(mps.empty());

  const size_t left_boundary = FindLeftBoundary(mps);
  const size_t initial_center = left_boundary + 1;
  for (size_t i = 0; i < initial_center; i++) {
    mps.dealloc(i);
  }

  assert(sites_set[0].size() == 2);
  const size_t total_event_size = sites_set.size();
  assert(world.size() >= Ly);
  const size_t event_size_every_group = total_event_size / Ly;
  const size_t group = world.rank();
  MeasuRes<TenElemT> measure_res;
  if (group < Ly) {
    const size_t site1 = sites_set[group * event_size_every_group][0];
    std::vector<size_t> site2_set;
    if (group != 0) {
      site2_set.reserve(event_size_every_group);
    } else {
      site2_set.reserve(total_event_size);
    }
    for (size_t i = 0; i < event_size_every_group; i++) {
      site2_set.push_back(sites_set[group * event_size_every_group + i][1]);
    }
    measure_res = MeasureTwoSiteOpGroup(mps, initial_center, phys_ops1, phys_ops2, site1, site2_set);
  }
  if (group >= Ly) {
    std::cout << "warning: processor " << world.rank() << " are idle." << std::endl;
  }

  /*
  std::string res_file_basename_group = res_file_basename + std::to_string(group);
  if(group < Ly){
    DumpMeasuRes(measure_res, res_file_basename_group); //tmp code to avoid communication bug
  }*/

  if (group == 0) {
    for (size_t recv_group = 1; recv_group < Ly; recv_group++) {
      const size_t site1 = sites_set[recv_group * event_size_every_group][0];
      std::vector<TenElemT> recved_avgs;
      world.recv(recv_group, recv_group, recved_avgs);
      for (size_t i = 0; i < event_size_every_group; i++) {
        const size_t site2 = sites_set[recv_group * event_size_every_group + i][1];
        measure_res.push_back(MeasuResElem<TenElemT>({site1, site2}, recved_avgs[i]));
      }
    }
  } else if (group < Ly) {
    std::vector<TenElemT> avgs;
    avgs.reserve(event_size_every_group);
    for (size_t i = 0; i < event_size_every_group; i++) {
      avgs.push_back(measure_res[i].avg);
    }
    world.send(0, group, avgs);
  }
  if (group == 0) {
    DumpMeasuRes(measure_res, res_file_basename);
  }
  return measure_res;
}

/**
 * This function used to measure two point function in such case:
 *   1. no insertion operator
 *   2. the site of first operator is fixed
 *   3. memory need saved
 * For more details please see below parameter interpretations.
 * @tparam TenElemT
 * @tparam QNT
 * @param mps   only has data of site initial_center, which is consistent with disk data
 * @param initial_center the mps center in disk data
 * @param phys_ops1
 * @param phys_ops2
 * @param site1   suppose site1 > initial_center
 * @param site2_set  suppose site2_set[end] > site2_set[end-1] > ....> site2_set[0] > site1
 * @return
 */
template<typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureTwoSiteOpGroup(
    FiniteMPS<TenElemT, QNT> &mps,
    const size_t initial_center,
    const QLTensor<TenElemT, QNT> &phys_ops1,
    const QLTensor<TenElemT, QNT> &phys_ops2,
    const size_t site1,
    const std::vector<size_t> &site2_set
) {
  std::string mps_path = kMpsPath;//usual case
  //move the center to site1
  for (size_t j = initial_center; j < site1; j++) {
    mps.LoadTen(j + 1, GenMPSTenName(mps_path, j + 1));
    mps.LeftCanonicalizeTen(j);
    mps.dealloc(j);
  }

  //Contract mps[site1]*phys_ops1*dag(mps[site1])
  auto id_op_set = mps.GetSitesInfo().id_ops;
  //Contract on site1
  std::vector<size_t> head_mps_ten_ctrct_axes1{1};
  std::vector<size_t> head_mps_ten_ctrct_axes2{0, 2};
  std::vector<size_t> head_mps_ten_ctrct_axes3{0, 1};
  QLTensor<TenElemT, QNT> temp_ten0;
  auto ptemp_ten = new QLTensor<TenElemT, QNT>;//TODO: delete
  Contract(
      &mps[site1], &phys_ops1,
      {{1}, {0}},
      &temp_ten0
  );
  QLTensor<TenElemT, QNT> mps_ten_dag = Dag(mps[site1]);
  Contract(
      &temp_ten0, &mps_ten_dag,
      {head_mps_ten_ctrct_axes2, head_mps_ten_ctrct_axes3},
      ptemp_ten
  );
  mps_ten_dag.GetBlkSparDataTen().Clear();//Save memory
  mps.dealloc(site1);

  size_t eated_site = site1; //the last site has been contracted
  MeasuRes<TenElemT> measure_res(site2_set.size());
  for (size_t event = 0; event < site2_set.size(); event++) {
    const size_t site2 = site2_set[event];
    while (eated_site < site2 - 1) {
      size_t eating_site = eated_site + 1;
      mps.LoadTen(eating_site, GenMPSTenName(mps_path, eating_site));
      //Contract ptemp_ten*mps[eating_site]*dag(mps[eating_site])
      CtrctMidTen(mps, eating_site, id_op_set[eating_site], id_op_set[eating_site], ptemp_ten);

      eated_site = eating_site;
      mps.dealloc(eated_site);
    }
    //now site2-1 has been eaten.
    mps.LoadTen(site2, GenMPSTenName(mps_path, site2));
    //Contract ptemp_ten*mps[site2]*ops2*dag(mps[site2]) gives the expected value.
    std::vector<size_t> tail_mps_ten_ctrct_axes1{0, 1, 2};
    std::vector<size_t> tail_mps_ten_ctrct_axes2{2, 0, 1};
    QLTensor<TenElemT, QNT> temp_ten2, temp_ten3, res_ten;
    Contract(&mps[site2], ptemp_ten, {{0}, {0}}, &temp_ten2);
    Contract(&temp_ten2, &phys_ops2, {{0}, {0}}, &temp_ten3);
    mps_ten_dag = Dag(mps[site2]);
    Contract(
        &temp_ten3, &mps_ten_dag,
        {tail_mps_ten_ctrct_axes1, tail_mps_ten_ctrct_axes2},
        &res_ten
    );
    measure_res[event] = MeasuResElem<TenElemT>({site1, site2}, res_ten());

    mps.dealloc(site2);//according now code this site2 will load again in next loop. This may be optimized one day.
  }
  delete ptemp_ten;
  return measure_res;
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
template<typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureElectronPhonon4PointFunction(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::vector<QLTensor<TenElemT, QNT>> &phys_ops,
    const std::vector<std::vector<size_t>> &sites_set,
    const std::string &res_file_basename
) {
  assert(mps.empty());
  const size_t left_boundary = FindLeftBoundary(mps);
  for (size_t i = 0; i <= left_boundary + 1; i++) {
    mps.dealloc(i);
  }
  const size_t Ly = 4; //change here when Ly of system changes
  std::cout << "note: Ly = " << Ly << std::endl;
  assert(phys_ops.size() == 4); // 4-point function
  assert(sites_set.size() % Ly == 0);
  for (size_t i = 0; i < sites_set.size(); i++) {
    assert(sites_set[i].size() == 4); // 4-point function
  }
  const size_t measure_event_num = sites_set.size();
  const size_t measure_event_per_group_num = sites_set.size() / Ly;
  for (size_t group = 0; group < Ly; group++) {
    std::vector<std::vector<size_t>> sites_set_group(
        sites_set.cbegin() + group * measure_event_per_group_num,
        sites_set.cbegin() + (group + 1) * measure_event_per_group_num
    );
    size_t site1 = sites_set_group[0][0];
    size_t site2 = sites_set_group[0][1];
    assert(sites_set_group.size() == measure_event_per_group_num);
    for (auto iter = sites_set_group.begin(); iter < sites_set_group.cend(); iter++) {
      assert((*iter)[0] == site1);
      assert((*iter)[1] == site2);
      //TODO: check order: a. site1<site2<3<4 b. for every 3,4, ascending order
    }
  }
  MeasuRes<TenElemT> measure_res;
  measure_res.reserve(measure_event_num);
  for (size_t group = 0; group < Ly; group++) {
    std::vector<std::vector<size_t>> sites_set_group(
        sites_set.cbegin() + group * measure_event_per_group_num,
        sites_set.cbegin() + (group + 1) * measure_event_per_group_num
    );
    auto measure_res_group = MeasureElectronPhonon4PointFunctionGroup(mps,
                                                                      phys_ops,
                                                                      sites_set_group,
                                                                      left_boundary + 1
    );
    measure_res.insert(measure_res.end(), measure_res_group.begin(), measure_res_group.end());
  }
  DumpMeasuRes(measure_res, res_file_basename);
  return measure_res;
}

template<typename TenElemT, typename QNT>
struct TempTensorWithContract2Ops {
  using TenT = QLTensor<TenElemT, QNT>;
  size_t idx; // The last site having be contract
  TenT *tmp_tensor;

  TempTensorWithContract2Ops(const size_t &idx, TenT *const &tmp_tensor) :
      idx(idx), tmp_tensor(tmp_tensor) {}

  void MoveOnTo(
      FiniteMPS<TenElemT, QNT> &mps,
      size_t site_to
  ) {
    assert(site_to >= idx);
    for (idx = idx + 1; idx <= site_to; idx++) {
      if (mps(idx) == nullptr) {
        mps.LoadTen(idx, GenMPSTenName(kMpsPath, idx));
      }
      CtrctMidTen(mps, idx, TenT(), TenT(), tmp_tensor);
      mps.dealloc(idx);
    }
    idx--;
  }
};

template<typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureElectronPhonon4PointFunctionGroup(
    FiniteMPS<TenElemT, QNT> &mps, //input and output mps is empty
    const std::vector<QLTensor<TenElemT, QNT>> &phys_ops,
    const std::vector<std::vector<size_t>> &sites_set,
    const size_t initial_center
) {
  assert(mps.empty());
  std::string mps_path = kMpsPath;
  const size_t site1 = sites_set[0][0];
  const size_t site2 = sites_set[0][1];
//  const size_t fermion_op_dim = mps[site1].GetIndexes()[1].dim();// should be 4 for SSH-Hubbard
//  const size_t bonson_op_dim = mps[site1+1].GetIndexes()[1].dim();// usually should be 2;
  const size_t fermion_op_dim(4), bonson_op_dim(2);
  assert(fermion_op_dim != bonson_op_dim);//we use dimension to differentiate boson site or fermion site
  using Tensor = QLTensor<TenElemT, QNT>;

  static bool is_f_initial = false;
  static Tensor f;
  if (!is_f_initial) {
    mps.LoadTen(0, GenMPSTenName(mps_path, 0));
    Index<QNT> index_out_fermion = mps[0].GetIndexes()[1];
    mps.dealloc(0);
    Index<QNT> index_in_fermion = InverseIndex(index_out_fermion);
    f = Tensor({index_in_fermion, index_out_fermion});
    f({0, 0}) = 1;
    f({1, 1}) = -1;
    f({2, 2}) = -1;
    f({3, 3}) = 1;
    is_f_initial = true;
  }
  mps.LoadTen(initial_center, GenMPSTenName(mps_path, initial_center));
  for (size_t j = initial_center; j < site1; j++) {
    mps.LoadTen(j + 1, GenMPSTenName(mps_path, j + 1));
    mps.LeftCanonicalizeTen(j);
    mps.dealloc(j);
  }

  auto id_op_set = mps.GetSitesInfo().id_ops;
  //Contract on site1
  std::vector<size_t> head_mps_ten_ctrct_axes1{1};
  std::vector<size_t> head_mps_ten_ctrct_axes2{0, 2};
  std::vector<size_t> head_mps_ten_ctrct_axes3{0, 1};
  QLTensor<TenElemT, QNT> temp_ten0;
  auto ptemp_ten = new QLTensor<TenElemT, QNT>;//delete when first called MoveOnTo
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

  mps.dealloc(site1);

  for (size_t i = site1 + 1; i < site2; ++i) {
    mps.LoadTen(i, GenMPSTenName(mps_path, i));
    if (mps[i].GetIndexes()[1].dim() == bonson_op_dim) {
      CtrctMidTen(mps, i, id_op_set[i], id_op_set[i], ptemp_ten);
    } else {
      CtrctMidTen(mps, i, f, id_op_set[i], ptemp_ten);
    }
    mps.dealloc(i);
  }
  mps.LoadTen(site2, GenMPSTenName(mps_path, site2));
  CtrctMidTen(mps, site2, phys_ops[1], id_op_set[site2], ptemp_ten);
  mps.dealloc(site2);
  TempTensorWithContract2Ops<TenElemT, QNT> tmp_tensor_with_c2o(site2, ptemp_ten);

  MeasuRes<TenElemT> measure_res(sites_set.size());
  for (size_t event = 0; event < sites_set.size(); event++) {//event means measure event
    size_t site3 = sites_set[event][2];
    size_t site4 = sites_set[event][3];
    tmp_tensor_with_c2o.MoveOnTo(mps, site3 - 1);//When first move on, ptemp_ten was deleted.
    ptemp_ten = new Tensor(*tmp_tensor_with_c2o.tmp_tensor);// deep copy
    mps.LoadTen(site3, GenMPSTenName(mps_path, site3));
    CtrctMidTen(mps, site3, phys_ops[2], id_op_set[site3], ptemp_ten);
    for (size_t i = site3 + 1; i < site4; ++i) {
      mps.LoadTen(i, GenMPSTenName(mps_path, i));
      if (mps[i].GetIndexes()[1].dim() == bonson_op_dim) {
        CtrctMidTen(mps, i, id_op_set[i], id_op_set[i], ptemp_ten);
      } else {
        CtrctMidTen(mps, i, f, id_op_set[i], ptemp_ten);
      }
    }
    // Deal with tail tensor.
    std::vector<size_t> tail_mps_ten_ctrct_axes1{0, 1, 2};
    std::vector<size_t> tail_mps_ten_ctrct_axes2{2, 0, 1};
    QLTensor<TenElemT, QNT> temp_ten2, temp_ten3, res_ten;
    mps.LoadTen(site4, GenMPSTenName(mps_path, site4));
    Contract(&mps[site4], ptemp_ten, {{0}, {0}}, &temp_ten2);
    delete ptemp_ten;
    Contract(&temp_ten2, &phys_ops.back(), {{0}, {0}}, &temp_ten3);
    mps_ten_dag = std::move(Dag(mps[site4]));
    Contract(
        &temp_ten3, &mps_ten_dag,
        {tail_mps_ten_ctrct_axes1, tail_mps_ten_ctrct_axes2},
        &res_ten
    );
    measure_res[event] = MeasuResElem<TenElemT>(sites_set[event], res_ten());
  }
  //clean the data
  size_t site3 = sites_set.back()[2];
  size_t site4 = sites_set.back()[3];
  for (size_t i = site3; i <= site4; i++) {
    mps.dealloc(i);
  }
  delete tmp_tensor_with_c2o.tmp_tensor;
  return measure_res;
}

}//qlmps

#endif // MY_MEASURE_H