#ifndef GQDOUBLE_H
#define GQDOUBLE_H
#include "gqten/gqten.h"

using gqten::QNCard;
using gqten::U1QNVal;
using gqten::GQTenIndexDirType;

using TenElemT = gqten::GQTEN_Double;
using U1U1QN = gqten::QN<U1QNVal, U1QNVal>;
using Tensor = gqten::GQTensor<TenElemT, U1U1QN>;

using QNSctT2 = gqten::QNSector<U1U1QN>;
using IndexT2 = gqten::Index<U1U1QN>;



const auto qn0 = U1U1QN(
	{QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}
    );


const IndexT2 pb_outF = IndexT2({
      QNSctT2(U1U1QN({QNCard("N", U1QNVal(2)), QNCard("Sz", U1QNVal( 0))}), 1),
      QNSctT2(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal( 1))}), 1),
      QNSctT2(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(-1))}), 1),
      QNSctT2(U1U1QN({QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal( 0))}), 1) },
      GQTenIndexDirType::OUT
    );
const auto pb_inF = gqten::InverseIndex(pb_outF);

const IndexT2 pb_outB = IndexT2({QNSctT2(qn0,2)}, GQTenIndexDirType::OUT);
const IndexT2 pb_inB = gqten::InverseIndex(pb_outB);

void OperatorInitial();

#endif