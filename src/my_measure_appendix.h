//
// Created by Hao-Xin on 2023/6/12.
//

#ifndef MY_MEASURE_APPENDIX_H
#define MY_MEASURE_APPENDIX_H

#include "my_measure.h"

namespace qlmps {
using namespace qlten;

/**
 * This function used to measure two point function in such case:
 *   1. with insertion operator
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
inline MeasuRes<TenElemT> MeasureTwoSiteFermionOpGroup(
    FiniteMPS<TenElemT, QNT> &mps,
    const size_t initial_center,
    const Qltensor<TenElemT, QNT> &phys_ops1,
    const Qltensor<TenElemT, QNT> &phys_ops2,
    const size_t site1,
    const std::vector<size_t> &site2_set
) {
  std::string mps_path = kMpsPath;//usual case
  const size_t bonson_op_dim(2);
  using Tensor = Qltensor<TenElemT, QNT>;

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
  Qltensor<TenElemT, QNT> temp_ten0;
  auto ptemp_ten = new Qltensor<TenElemT, QNT>;//TODO: delete
  Contract(
      &mps[site1], &phys_ops1,
      {{1}, {0}},
      &temp_ten0
  );
  Qltensor<TenElemT, QNT> mps_ten_dag = Dag(mps[site1]);
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
      if (mps[eating_site].GetIndexes()[1].dim() == bonson_op_dim) {
        CtrctMidTen(mps, eating_site, id_op_set[eating_site], id_op_set[eating_site], ptemp_ten);
      } else {
        CtrctMidTen(mps, eating_site, f, id_op_set[eating_site], ptemp_ten);
      }

      //Contract ptemp_ten*mps[eating_site]*dag(mps[eating_site])

      eated_site = eating_site;
      mps.dealloc(eated_site);
    }
    //now site2-1 has been eaten.
    mps.LoadTen(site2, GenMPSTenName(mps_path, site2));
    //Contract ptemp_ten*mps[site2]*ops2*dag(mps[site2]) gives the expected value.
    std::vector<size_t> tail_mps_ten_ctrct_axes1{0, 1, 2};
    std::vector<size_t> tail_mps_ten_ctrct_axes2{2, 0, 1};
    Qltensor<TenElemT, QNT> temp_ten2, temp_ten3, res_ten;
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
Measure a two-site operator with insertion operator.
MPI version, memory optimized version

@tparam TenElemT Type of the tensor element, real or complex.
@param mps To-be-measured MPS.
@param sites_set The indexes of the two physical operators with ascending order
       for each measure event. Its size defines the number of measure events.
@param res_file_basename The basename of the output file.
*/
template<typename TenElemT, typename QNT>
inline MeasuRes<TenElemT> MeasureTwoSiteFermionOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const Qltensor<TenElemT, QNT> &phys_ops1,
    const Qltensor<TenElemT, QNT> &phys_ops2,
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
    measure_res = MeasureTwoSiteFermionOpGroup(mps, initial_center, phys_ops1, phys_ops2, site1, site2_set);
  }
  if (group >= Ly) {
    std::cout << "warning: processor " << world.rank() << " are idle." << std::endl;
  }

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

template<typename TenElemT, typename QNT>
MeasuRes<TenElemT> MeasureOnePhoneOp(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::vector<Qltensor<TenElemT, QNT>> &op_vec,
    const std::vector<size_t> &boson_set, // pseudo-sites
    const std::string &res_file_basename
) {
  size_t N = mps.size();
  size_t Nps = op_vec.size(); // pseudo-site number
  size_t N_ph = boson_set.size() / Nps;
  MeasuRes<TenElemT> measu_res;
  measu_res.reserve(N_ph);

  const size_t left_boundary = FindLeftBoundary(mps);
  const size_t initial_center = left_boundary + 1;

  for (size_t i = initial_center; i > 1; i--) { // canonicalize to first pseudo-site of phonon
    mps.RightCanonicalizeTen(i);
  }

  const std::string mps_path = kMpsPath;
  size_t last_boson_site_start = 0;
  for (size_t i = 0; i < N_ph; i++) {
    std::cout << "phonon site = " << i << std::endl;
    size_t boson_site_start = boson_set[i * Nps];
    assert(boson_site_start < N);
    if (i == 0) {
      for (size_t site = boson_site_start; site < boson_site_start + Nps; site++) {
        if (mps(site) == nullptr) {
          mps.LoadTen(site, GenMPSTenName(mps_path, site));
        }
      }
    } else {
      for (size_t site = last_boson_site_start; site < boson_site_start + Nps; site++) {
        if (mps(site) == nullptr) {
          mps.LoadTen(site, GenMPSTenName(mps_path, site));
        }
      }
      for (size_t site = last_boson_site_start; site < boson_site_start; site++) {
        mps.LeftCanonicalizeTen(site);
        mps.dealloc(site);
      }
    }
    TenElemT avg = OpsVecAvg(mps, op_vec, boson_site_start, boson_site_start + Nps - 1);
    measu_res.push_back(MeasuResElem<TenElemT>({i}, avg));
    last_boson_site_start = boson_site_start;
  }
  mps.clear();
  DumpMeasuRes(measu_res, res_file_basename);
  return measu_res;
}

}
#endif //MY_MEASURE_APPENDIX_H
