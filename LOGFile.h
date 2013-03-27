#ifndef INCLUDE_LOGFILE
#define INCLUDE_LOGFILE

#include "IOFile.h"

/*
1. User can define specified class derived from IOFile.

2. In the source code, the line number of txtfile start from '1', which is consistent with vi. 
   Be caution, the index of array starts from '0'.

3. Be very careful, when you use operator '>>' to read the final word of a line, 
   tellg() will give the position of last non-whitespace char, not the whitespace char (eg. 'endl').
   So the file read pointer still at the current line, not the next line.

4. On the oppsite, getline() will read the whole line including the 'endl' char,
   and then the file read pointer will be at next line now.

5. Last updated 2010-02-16 by YF@LiGroup, UWChem.
*/


// *** class CLOGFileReadIn *** //
class CLOGFileReadIn: public IOFile
{
public:
  bool ReadInAtoms();
  bool ReadInBasisAO();
  bool ReadInOverlap();
  bool ReadInMOCoeffAlph();
  bool ReadInMOCoeffBeta();
  
  void PrintAtoms();
  void PrintBasisAO();
  void PrintOverlap();
  void PrintMOCoeffAlph();
  void PrintMOCoeffBeta();
  void PrintMOEigenAlph();
  void PrintMOEigenBeta();
  
  int ReturnNAtoms(){ return nAtoms;   };
  int ReturnNBasis(){ return nAOBasis; };
  int ReturnNElecA(){ return nElecA;   };
  int ReturnNElecB(){ return nElecB;   };
  
  const double ( *(ReturnOverlap)( )    )[nBASIS_MAX]{ return overlapMatrix; };    //return 2-D array 
  const double ( *(ReturnMOCoeffAlph)() )[nBASIS_MAX]{ return coeffMOAlph;   };
  const double ( *(ReturnMOCoeffBeta)() )[nBASIS_MAX]{ return coeffMOBeta;   };
  
  const double* ReturnMOEigenAlph(){ return eigenAlph; };                          //return 1-D array
  const double* ReturnMOEigenBeta(){ return eigenBeta; };
  
  const CAtoms*      ReturnAtoms(){ return atoms;   };
  const CBasisSetAO* ReturnAOBasis(){ return AOBasis; };
  
private:
  void ReadInMOCoeff( double (*coeff)[nBASIS_MAX], double* eigen  );
  
  int nAtoms;
  CAtoms atoms[nATOMS_MAX];
  
  int nAOBasis;
  CBasisSetAO AOBasis[nBASIS_MAX];
  
  int nElecA, nElecB;
  int nVirtA, nVirtB;
  
  double overlapMatrix[nBASIS_MAX][nBASIS_MAX];      //Overlap matrix
  
  double coeffMOAlph[nBASIS_MAX][nBASIS_MAX];        //MO coeff matrix
  double coeffMOBeta[nBASIS_MAX][nBASIS_MAX];
  
  double eigenAlph[nBASIS_MAX];                      //Eigenvalue marix, au
  double eigenBeta[nBASIS_MAX];
};

#endif

