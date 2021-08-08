#include "gqdouble.h"

Tensor sz = Tensor({pb_inF, pb_outF});
Tensor sp = Tensor({pb_inF, pb_outF});
Tensor sm = Tensor({pb_inF, pb_outF});
Tensor id = Tensor({pb_inF, pb_outF});
Tensor sx = Tensor({pb_inF, pb_outF});
Tensor isy = Tensor({pb_inF, pb_outF}); // 1i*sy to write as a real value, note to change the coefficient,i.e. sy*sy=- isy*isy
//auto sy = Tensor({pb_inF, pb_outF});

Tensor f = Tensor({pb_inF,pb_outF}); //fermion's insertion operator

auto bupc = Tensor({pb_inF,pb_outF}); //hardcore boson, b_up^creation, used for JW transformation
auto bupa = Tensor({pb_inF,pb_outF}); //hardcore boson, b_up^annihilation
auto bdnc = Tensor({pb_inF,pb_outF}); //hardcore boson, b_down^creation
auto bdna = Tensor({pb_inF,pb_outF}); //hardcore boson, b_down^annihilation


auto bupcF = Tensor({pb_inF,pb_outF}); // matrix product of bupc * f
auto bupaF = Tensor({pb_inF,pb_outF});
auto Fbdnc = Tensor({pb_inF,pb_outF});
auto Fbdna = Tensor({pb_inF,pb_outF});


auto cupccdnc = Tensor({pb_inF,pb_outF}); // c_up^creation * c_down^creation=b_up^creation*b_down^creation*F

auto cdnacupa = Tensor({pb_inF,pb_outF}); // onsite pair, usually c_up*c_dn


auto Uterm = Tensor({pb_inF,pb_outF}); // Hubbard Uterm, nup*ndown

auto nf = Tensor({pb_inF,pb_outF}); // nup+ndown, fermion number

auto nfsquare = Tensor({pb_inF,pb_outF}); // nf^2
auto nup = Tensor({pb_inF,pb_outF}); // fermion number of spin up
auto ndn = Tensor({pb_inF,pb_outF}); // ndown


// Bosonic operators
auto idB = Tensor({pb_inB, pb_outB}); //

auto a_create = Tensor({pb_inB, pb_outB}); // phonon creation operator
auto a_annihilate = Tensor({pb_inB, pb_outB});

auto &a =a_annihilate;
auto &adag = a_create;
auto n_a =  Tensor({pb_inB, pb_outB}); // the number of phonon

auto x = Tensor({pb_inB, pb_outB}); // a_creation + a_annihilate
auto aadag = Tensor({pb_inB, pb_outB});
void OperatorInitial(){
    static bool initialized = false;
    if(!initialized){
        sz({1, 1}) = 0.5;
        sz({2, 2}) = -0.5;
        sp({1, 2}) = 1.0;
        sm({2, 1}) = 1.0;
        id({0,0}) = 1;
        id({1,1}) = 1;
        id({2,2}) = 1;
        id({3,3}) = 1;
        sx({1,2}) = 0.5;
        sx({2,1}) = 0.5;
        isy({1,2}) = 0.5;
        isy({2,1}) = -0.5;
        //sy({1,2}) = -0.5i;
        //sy({2,1}) = 0.5i;

        f({0,0}) = 1;
        f({1,1}) = -1;
        f({2,2}) = -1;
        f({3,3}) = 1;


        bupc({0,2}) = 1;
        bupc({1,3}) = 1;
        bdnc({0,1}) = 1;
        bdnc({2,3}) = 1;
        bupa({2,0}) = 1;
        bupa({3,1}) = 1;
        bdna({1,0}) = 1;
        bdna({3,2}) = 1;


        bupcF({0,2}) = -1;
        bupcF({1,3}) = 1;
        Fbdnc({0,1}) = 1;
        Fbdnc({2,3}) = -1;
        bupaF({2,0}) = 1;
        bupaF({3,1}) = -1;
        Fbdna({1,0}) = -1;
        Fbdna({3,2}) = 1;


        cupccdnc({0,3}) = 1;
        cdnacupa({3,0}) = 1;


        Uterm({0,0}) = 1;


        nf({0,0}) = 2;
        nf({1,1}) = 1;
        nf({2,2}) = 1;


        nfsquare({0,0}) = 4;
        nfsquare({1,1}) = 1;
        nfsquare({2,2}) = 1;

        nup({0,0}) = 1;
        nup({1,1}) = 1;
        ndn({0,0}) = 1;
        ndn({2,2}) = 1;


        idB({0,0}) = 1;
        idB({1,1}) = 1;

        a_create({0,1}) = 1;
        a_annihilate({1,0}) = 1;

        n_a({0,0}) = 1;

        x({0,1}) = 1;
        x({1,0}) = 1;
        aadag = idB+(-n_a);
        initialized=true;
    }
}
