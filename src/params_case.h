#pragma once
#include "gqmps2/gqmps2.h"
using gqmps2::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Lx = ParseInt("Lx");
    Ly = ParseInt("Ly");
    t = ParseDouble("t");
    g = ParseDouble("g");
    U = ParseDouble("U");
    Np = ParseInt("Np");
    noise = ParseDoubleVec("noise");
    omega = ParseDouble("omega");
    Numhole = ParseInt("Numhole");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    TotalThreads = ParseInt("TotalThreads");
    SvdOuterThreads = ParseInt("SvdOuterThreads");
    Perturbation = ParseBool("Perturbation");
    if (Perturbation) {
      PA = ParseDouble("PerturbationAmplitude");
      PerturbationPeriod = ParseInt("PerturbationPeriod");
    } else {
      PA = 0.0;
      PerturbationPeriod = 1;
    }

  }

  long Lx;
  long Ly;
  long Np; // The number of pseudosite
  double t;
  double g;
  double U;
  double omega;
  long Sweeps;
  long Dmin;
  long Dmax;
  int Numhole;
  double CutOff;
  double LanczErr;
  long MaxLanczIter;
  int TotalThreads;
  int SvdOuterThreads;
  std::vector<double> noise;
  bool Perturbation;
  double PA;
  size_t PerturbationPeriod;
};