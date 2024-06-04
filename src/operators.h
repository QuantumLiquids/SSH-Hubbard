

#ifndef PHONON_ELECTRON_OPERATORS2_H
#define PHONON_ELECTRON_OPERATORS2_H
#include "qlten/qlten.h"

using TenElemT = qlten::QLTEN_Double;
using U1U1QN = qlten::QN<U1QNVal, U1QNVal>;
using Tensor = qlten::Qltensor<TenElemT, U1U1QN>;


//Fermionic operators
extern Tensor sz, sp,sm, id ,sx, isy;//sy
extern Tensor f, bupc, bupa, bdnc, bdna;
extern Tensor bupcF, bupaF, Fbdnc, Fbdna;
extern Tensor cupccdnc, cdnacupa, Uterm, nf, nfsquare, nup,ndn;

// Bosonic operators
extern Tensor idB, a_create, a_annihilate;
extern Tensor &a ;//=a_annihilate;
extern Tensor &adag; //= a_create;
extern Tensor n_a, aadag, x;
#endif