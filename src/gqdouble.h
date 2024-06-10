#ifndef GQDOUBLE_H
#define GQDOUBLE_H
#include "qlten/qlten.h"

using qlten::QNCard;
using qlten::U1QNVal;
using qlten::TenIndexDirType;

using TenElemT = qlten::QLTEN_Double;
using U1U1QN = qlten::QN<U1QNVal, U1QNVal>;
using Tensor = qlten::QLTensor<TenElemT, U1U1QN>;

using QNSctT = qlten::QNSector<U1U1QN>;
using IndexT = qlten::Index<U1U1QN>;

const auto qn0 = U1U1QN(
    {QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}
);

const IndexT pb_outF = IndexT({QNSctT(U1U1QN({QNCard("N", U1QNVal(2)), QNCard("Sz", U1QNVal(0))}), 1),
                               QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(1))}), 1),
                               QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(-1))}), 1),
                               QNSctT(U1U1QN({QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}), 1)},
                              TenIndexDirType::OUT
);
const auto pb_inF = qlten::InverseIndex(pb_outF);

const IndexT pb_outB = IndexT({QNSctT(qn0, 2)}, TenIndexDirType::OUT);
const IndexT pb_inB = qlten::InverseIndex(pb_outB);

void OperatorInitial();

#endif