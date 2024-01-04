#ifndef GQDOUBLE_H
#define GQDOUBLE_H
#include "gqten/gqten.h"

using gqten::QNCard;
using gqten::U1QNVal;
using gqten::GQTenIndexDirType;

using TenElemT = gqten::GQTEN_Double;
using U1U1QN = gqten::QN<U1QNVal, U1QNVal>;
using Tensor = gqten::GQTensor<TenElemT, U1U1QN>;

using QNSctT = gqten::QNSector<U1U1QN>;
using IndexT = gqten::Index<U1U1QN>;

const auto qn0 = U1U1QN(
    {QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}
);

const IndexT pb_outF = IndexT({QNSctT(U1U1QN({QNCard("N", U1QNVal(2)), QNCard("Sz", U1QNVal(0))}), 1),
                               QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(1))}), 1),
                               QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(-1))}), 1),
                               QNSctT(U1U1QN({QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}), 1)},
                              GQTenIndexDirType::OUT
);
const auto pb_inF = gqten::InverseIndex(pb_outF);

const IndexT pb_outB = IndexT({QNSctT(qn0, 2)}, GQTenIndexDirType::OUT);
const IndexT pb_inB = gqten::InverseIndex(pb_outB);

void OperatorInitial();

#endif