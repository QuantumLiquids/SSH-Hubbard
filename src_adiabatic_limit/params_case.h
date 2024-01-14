#pragma once
#include "gqmps2/gqmps2.h"
using gqmps2::CaseParamsParserBasic;

enum PhononConfigIntialType {
  Zero = 0,
  Random
};

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    t = ParseDouble("t");
    U = ParseDouble("U");
    alpha = ParseDouble("alpha");
    noise = ParseDoubleVec("noise");
    K = ParseDouble("K");
    Numhole = ParseInt("Numhole");
    HamiltonianIters = ParseInt("HamiltonianIters");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    TotalThreads = ParseInt("TotalThreads");
    Perturbation = ParseBool("Perturbation");
    if (Perturbation) {
      PA = ParseDouble("PerturbationAmplitude");
      PerturbationPeriod = ParseInt("PerturbationPeriod");
    } else {
      PA = 0.0;
      PerturbationPeriod = 1;
    }
    ph_init_config = static_cast<PhononConfigIntialType>(ParseInt("InitialPhononConfig"));
  }

  size_t Lx;
  size_t Ly;
  double t;
  double U;
  double alpha;
  double K;
  size_t HamiltonianIters;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  int Numhole;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  size_t TotalThreads;
  PhononConfigIntialType ph_init_config;
  std::vector<double> noise;
  bool Perturbation;
  double PA;
  size_t PerturbationPeriod;
};