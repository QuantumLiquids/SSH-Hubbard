#include "gqdouble.h"

Tensor sz = Tensor({pb_in, pb_out});
Tensor sp = Tensor({pb_in, pb_out});
Tensor sm = Tensor({pb_in, pb_out});
Tensor id = Tensor({pb_in, pb_out});
Tensor sx = Tensor({pb_in, pb_out});
Tensor isy =
    Tensor({pb_in, pb_out}); // 1i*sy to write as a real value, note to change the coefficient,i.e. sy*sy=- isy*isy
//auto sy = Tensor({pb_in, pb_out});

Tensor f = Tensor({pb_in, pb_out}); //fermion's insertion operator

auto bupc = Tensor({pb_in, pb_out}); //hardcore boson, b_up^creation, used for JW transformation
auto bupa = Tensor({pb_in, pb_out}); //hardcore boson, b_up^annihilation
auto bdnc = Tensor({pb_in, pb_out}); //hardcore boson, b_down^creation
auto bdna = Tensor({pb_in, pb_out}); //hardcore boson, b_down^annihilation


auto bupcF = Tensor({pb_in, pb_out}); // matrix product of bupc * f
auto bupaF = Tensor({pb_in, pb_out});
auto Fbdnc = Tensor({pb_in, pb_out});
auto Fbdna = Tensor({pb_in, pb_out});

auto cupccdnc = Tensor({pb_in, pb_out}); // c_up^creation * c_down^creation=b_up^creation*b_down^creation*F

auto cdnacupa = Tensor({pb_in, pb_out}); // onsite pair, usually c_up*c_dn


auto Uterm = Tensor({pb_in, pb_out}); // Hubbard Uterm, nup*ndown

auto nf = Tensor({pb_in, pb_out}); // nup+ndown, fermion number

auto nfsquare = Tensor({pb_in, pb_out}); // nf^2
auto nup = Tensor({pb_in, pb_out}); // fermion number of spin up
auto ndn = Tensor({pb_in, pb_out}); // ndown


void OperatorInitial() {
  static bool initialized = false;
  if (!initialized) {
    sz({1, 1}) = 0.5;
    sz({2, 2}) = -0.5;
    sp({1, 2}) = 1.0;
    sm({2, 1}) = 1.0;
    id({0, 0}) = 1;
    id({1, 1}) = 1;
    id({2, 2}) = 1;
    id({3, 3}) = 1;
    sx({1, 2}) = 0.5;
    sx({2, 1}) = 0.5;
    isy({1, 2}) = 0.5;
    isy({2, 1}) = -0.5;
    //sy({1,2}) = -0.5i;
    //sy({2,1}) = 0.5i;

    f({0, 0}) = 1;
    f({1, 1}) = -1;
    f({2, 2}) = -1;
    f({3, 3}) = 1;

    bupc({0, 2}) = 1;
    bupc({1, 3}) = 1;
    bdnc({0, 1}) = 1;
    bdnc({2, 3}) = 1;
    bupa({2, 0}) = 1;
    bupa({3, 1}) = 1;
    bdna({1, 0}) = 1;
    bdna({3, 2}) = 1;

    bupcF({0, 2}) = -1;
    bupcF({1, 3}) = 1;
    Fbdnc({0, 1}) = 1;
    Fbdnc({2, 3}) = -1;
    bupaF({2, 0}) = 1;
    bupaF({3, 1}) = -1;
    Fbdna({1, 0}) = -1;
    Fbdna({3, 2}) = 1;

    cupccdnc({0, 3}) = 1;
    cdnacupa({3, 0}) = 1;

    Uterm({0, 0}) = 1;

    nf({0, 0}) = 2;
    nf({1, 1}) = 1;
    nf({2, 2}) = 1;

    nfsquare({0, 0}) = 4;
    nfsquare({1, 1}) = 1;
    nfsquare({2, 2}) = 1;

    nup({0, 0}) = 1;
    nup({1, 1}) = 1;
    ndn({0, 0}) = 1;
    ndn({2, 2}) = 1;

    initialized = true;
  }
}
