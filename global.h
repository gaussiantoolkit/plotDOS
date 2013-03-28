#ifndef INCLUDE_GLOBAL
#define INCLUDE_GLOBAL

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <limits>
#include <cfloat>
#include <cmath>
#include <ctime>

using namespace std;

//const
const int nATOMS_MAX=1000;
const int nBASIS_MAX=5000;
const int nTYPES_MAX= 100;

const double au2ev=27.2114;

//
struct CAtoms
{
  int center;
  int atomic;
  int atomic_type;
  double xyz[3];
};

//
struct CBasisSetAO
{
  int center;
  int atomic;
  int l;
  int m;
};

//
struct CCategoryAO
{
  int center;               //only meaningful for 'individual TM D type AO'
  int atomic;
  int l;
  int nAOs;
  int AO[nBASIS_MAX];       //the list of AO with the same type
};

//
const string str_L[] = { "s", "p", "d", "f" };
const string str_Atom[] = { " X",
                            " H","He","Li","Be"," B"," C"," N"," O"," F","Ne",
                            "Na","Mg","Al","Si"," P"," S","Cl","Ar"," K","Ca",
                            "Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                            "Ga","Ge","As","Se","Br","Kr","Rb","Sr"," Y","Zr",
                            "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
                            "Sb","Te"," I","Xe","Cs","Ba",
                          };

#endif

