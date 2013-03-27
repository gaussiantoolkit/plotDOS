#ifndef INCLUDE_DOS
#define INCLUDE_DOS

#include "global.h"

// ***  class CCataAOs  *** //
class CCataAOs
{
public:
  CCataAOs( int nAOs=0 ){ nBasis = nAOs; };
  
  void CategorizeAO( const CBasisSetAO* AO );    //categorize basisset according to AO Atomic types.
  void CategorizeTM( const CBasisSetAO* AO );    //only categorize TM's D type AOs.
  void PrintHead();
  void PrintAOCategory();
  
  bool IsTM(int i){ if(i>=21&&i<=29) return true; else return false; };
  
  int ReturnNAOTypes(){ return nAOTypes; };
  const CCategoryAO* ReturnAOTypes(){ return AOTypes; };
  
private:
  int FindAOType( int iAO, const CBasisSetAO* AO );
  int FindTMType( int iAO, const CBasisSetAO* AO );
  int nBasis;
  int nAOTypes;
  
  CCategoryAO AOTypes[nTYPES_MAX];
};


// *** class CCalcDOS *** //
class CCalcDOS
{
public:
  CCalcDOS( int nAOs=0, int nTypes=0 ){  nBasis = nAOs; nAOTypes = nTypes; };   //constructer
  
  bool CalcAOperMO( const double (*overlap)[nBASIS_MAX], const double(*coeff)[nBASIS_MAX], int iMO );    // calculate iMO
  void CalcAOperMO( const double (*overlap)[nBASIS_MAX], const double(*coeff)[nBASIS_MAX] );             // calculate all MO
  void PrintAOperMO();
  void CalcNormalAO();                                  //sum AO percent for each MOs.
  
  bool SumCTperMO( const CCategoryAO*, int iMO );       //sum AO percent according to AO types for 'iMO'
  void SumCTperMO( const CCategoryAO* );                //... for each MOs.
  void CalcNormalCT();                                  //sum CT percent for each MOs.
  
  void SumElecDensity( int nElec );                  //sum elec density
  void SumElecDensity( int iStart, int iEnd );       //sum elec density from iStart to iEnd, according to AO types.
  void PrintElecDensity( const CCategoryAO* );
  
  void PrintCTperMO( const CCategoryAO* );
  void PrintCTperMO( const CCategoryAO*, int iStart, int iEnd, int iHOMO = 0 );
  
  void SumDOSPartial(const double* eigen, double binStart = -20.0, double binEnd = 20.0 , double binSize = 0.2 );
  void PrintDOSPartial( const CCategoryAO* );
//KB start
  void PrintDOSPartialBeta( const CCategoryAO* );
//KB end
  void SumDOSTotal(const double* eigen, double binStart = -20.0, double binEnd = 20.0 , double binSize = 0.2 );
  void PrintDOSTotal( const CCategoryAO* );
//KB start
  void PrintDOSTotalBeta( const CCategoryAO* );
//KB end
  
  void FillEigen(const double *eigen);
private:
   
  int  FillEBins(double binStart, double binEnd, double binSize);
  
  int nBasis;
  int nAOTypes;
  
  double normalAO[nBASIS_MAX];                   //normalize const of each MO.      ="1"
  double normalCT[nTYPES_MAX];                   //normalize const of each type MO. ="n"
  
  double AOperMO[nBASIS_MAX][nBASIS_MAX];        //AO percentage of MO according to mulliken pop analysis.
  double CTperMO[nBASIS_MAX][nTYPES_MAX];        //AO percentage of MO,according to category.
  
  double eigenEV[nBASIS_MAX];                    //au 2 ev
  vector<double> eigenBins;
  vector<int>    eigenNDOS;                      //total MOs, == eigenPDOS[i]
  
  vector<double> eigenPDOS;                      //sum eigenCT[i][0-nAOTypes];
  vector< vector<double> >eigenCT;
  
  vector<int> eigenMOStart;                      //iMO Start
  vector<int> eigenMOEnd;                        //iMO End
  
  vector<double> elecDensity;                    //elec density
  
  
};

#endif

