#ifndef GQDOUBLE_H
#define GQDOUBLE_H
#include "gqten/gqten.h"

using gqten::QNCard;
using gqten::U1QNVal;
using gqten::GQTenIndexDirType;

using TenElemT = gqten::GQTEN_Double;
using U1U1QN = gqten::special_qn::U1U1QN;
using Tensor = gqten::GQTensor<TenElemT, U1U1QN>;

using QNSctT = gqten::QNSector<U1U1QN>;
using IndexT = gqten::Index<U1U1QN>;

const auto qn0 = U1U1QN(
    {QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}
);

const IndexT pb_out = IndexT({QNSctT(U1U1QN({QNCard("N", U1QNVal(2)), QNCard("Sz", U1QNVal(0))}), 1),
                              QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(1))}), 1),
                              QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(-1))}), 1),
                              QNSctT(U1U1QN({QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}), 1)},
                             GQTenIndexDirType::OUT
);
const auto pb_in = gqten::InverseIndex(pb_out);

void OperatorInitial();

#endif