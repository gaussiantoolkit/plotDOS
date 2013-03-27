#include "DOS.h"

//To test the AO[iAO] is existed at categorizeAO[i] (return i ) or not (return -1);
int CCataAOs::FindAOType( int iAO, const CBasisSetAO* AO )
{
  for( int i=0; i<nAOTypes; i++)
  {
    if( AO[iAO].atomic==AOTypes[i].atomic && AO[iAO].l==AOTypes[i].l )  return i;
  };
  
  return -1;   //this is the new AO type
};

//
int CCataAOs::FindTMType( int iAO, const CBasisSetAO* AO )
{
  for( int i=0; i<nAOTypes; i++)
  {
    if(   AO[iAO].atomic==AOTypes[i].atomic 
       && AO[iAO].l     ==AOTypes[i].l 
       && AO[iAO].center==AOTypes[i].center )
      return i;
  };
  
  return -1;   //this is the new AO type
};


//categorize AO
void CCataAOs::CategorizeAO( const CBasisSetAO* AO)
{
  int i,j,k;
  nAOTypes = 0;                                     //initialize
  for(i=0;i<nTYPES_MAX;i++)   
  {
    AOTypes[i].nAOs   =  0;
    AOTypes[i].center = -1;
  };
  
  //categorize AOs
  for(j=0;j<nBasis;j++)
  {
    k = FindAOType(j, AO);
    if( k == -1 )
    {
      //add a new AO category
      AOTypes[nAOTypes].atomic = AO[j].atomic;     //trick
      AOTypes[nAOTypes].l      = AO[j].l;
      AOTypes[nAOTypes].nAOs   = 1;
      AOTypes[nAOTypes].AO[0]  = j;
      nAOTypes++;
    }else
    {
      //update a exsiting AO category
      AOTypes[k].AO[AOTypes[k].nAOs] = j;         //trick
      AOTypes[k].nAOs++;
    };
  };
};

//categorize TM AO
void CCataAOs::CategorizeTM( const CBasisSetAO* AO)
{
  int i,j,k;
  nAOTypes = 0;                                      //initialize
  for(i=0;i<nTYPES_MAX;i++)
  {
    AOTypes[i].nAOs   =  0;
    AOTypes[i].center = -1;
  };
  
  //categorize AOs
  for(j=0;j<nBasis;j++)
  {
    if( IsTM(AO[j].atomic) && AO[j].l == 2 )
    {
      k = FindTMType(j, AO);
      if( k == -1 )
      {
        //add a new AO category
        AOTypes[nAOTypes].atomic = AO[j].atomic;     //trick
        AOTypes[nAOTypes].center = AO[j].center;
        AOTypes[nAOTypes].l      = AO[j].l;
        AOTypes[nAOTypes].nAOs   = 1;
        AOTypes[nAOTypes].AO[0]  = j;
        nAOTypes++;
      }else
      {
        //update a exsiting AO category
        AOTypes[k].AO[AOTypes[k].nAOs] = j;          //trick
        AOTypes[k].nAOs++;
      };
    };
  };
};

//print
void CCataAOs::PrintHead()
{
};

void CCataAOs::PrintAOCategory()
{
  cout << "nAO Types:" << nAOTypes << endl;
  cout << "--iType----atomic----L----nAOs--" << endl;
  for(int i=0;i<nAOTypes;i++)
  {
    cout << setw(6) << i+1  << " ";
    cout << setw(8) << AOTypes[i].atomic << " ";
    cout << setw(6) << AOTypes[i].l      << " ";
    cout << setw(6) << AOTypes[i].nAOs   << endl;
  };
  cout << "--------------------------------" << endl;
};


//calculate every AO contribution of 'iMO', based on mulikin pop.
bool CCalcDOS::CalcAOperMO( const double (*overlap)[nBASIS_MAX], const double(*coeff)[nBASIS_MAX], int iMO )
{
  if( iMO < 0 || iMO >= nBasis ) return false;
  int j,k;  
  for(j=0;j<nBasis;j++)
  {
    
    AOperMO[iMO][j] = 0;
    
    //mulikin pop with overlap matrix;
    //phi=c1*AO1+c2*AO2, the AO1 percentage of phi = C1*C1 + C1*C2*<AO1|AO2>
    for(k=0;k<nBasis;k++)  AOperMO[iMO][j] += coeff[iMO][j]*coeff[iMO][k]*overlap[j][k];
    
    //only count on the diagonal coeff. C1*C1
    //AOperMO[iMO][j] = coeff[iMO][j] * coeff[iMO][j];
  
  };
  return true;
};

//calculate all MO
void CCalcDOS::CalcAOperMO( const double (*overlap)[nBASIS_MAX], const double(*coeff)[nBASIS_MAX] )
{
  //jwm cout <<"Calculate AO per MO :: ";
  int i,j;
  int iDotStep = nBasis / 20 + 1; 
  
  for(i=0;i<nBasis;i++)
  { 
    if( i % iDotStep == 0 ); //jwm cout << ". " << flush;
    CalcAOperMO(overlap,coeff,i);                        //calling CalcAOperMO(), very time consuming
  };
  //jwm cout << endl;
  
};

//
void CCalcDOS::CalcNormalAO()
{
  for(int i=0;i<nBasis;i++)
  {
    normalAO[i] = 0;
    for(int j=0;j<nBasis;j++) normalAO[i] += AOperMO[i][j];
  };
};

//
void CCalcDOS::PrintAOperMO()
{
  int i,j;
  cout.precision(5);
  cout << "---iMO--Normalize--------Coeff--------" << endl;
  
  //print AOperMO[][]
  for(i=0;i<nBasis;i++)
  {
    cout << setw(5) << i+1 << " ";
    cout << setw(10) << fixed << normalAO[i] + 0.00001 << "     ";
    if( nBasis < 12 )
    {
      for(j=0;j<nBasis;j++)         cout << setw(8) << fixed << AOperMO[i][j] << " ";
    } else
    {
      for(j=0;       j<6;     j++)  cout << setw(8) << fixed << AOperMO[i][j] << " ";
      cout << " ...... ";
      for(j=nBasis-6;j<nBasis;j++)  cout << setw(8) << fixed << AOperMO[i][j] << " ";
    };
    
    cout << endl;
  };
  cout << "----------------------" << endl;
};

bool CCalcDOS::SumCTperMO( const CCategoryAO* AOCategory, int iMO )
{
  int j,k;
  for(j=0;j<nAOTypes;j++)
  {
    CTperMO[iMO][j] = 0;   //init
    for(k=0;k<AOCategory[j].nAOs;k++)   CTperMO[iMO][j] += AOperMO[iMO][AOCategory[j].AO[k]];
  };
};

//
void CCalcDOS::SumCTperMO( const CCategoryAO* AOCategory )
{
  for(int i=0;i<nBasis;i++)
    SumCTperMO( AOCategory, i );      //calling SumCTperMO();
};

//
void CCalcDOS::CalcNormalCT()
{
  for(int i=0;i<nBasis;i++)
  {
    normalCT[i] = 0;      //init
    for(int j=0;j<nAOTypes;j++) normalCT[i] += CTperMO[i][j];
  };
};

//
void CCalcDOS::PrintCTperMO( const CCategoryAO* AOCategory, int iStart, int iEnd, int iHOMO )
{
  int i,j;
  cout.precision(3);
  
  //print head
  cout << setw(5)<< "iMO" ;
  cout << setw(8)<< "E(ev)";
  
  for(i=0;i<nAOTypes;i++) 
  {
    if( AOCategory[i].center != -1 ) 
      // 'AOCategory[i].center' will be the label of this AO atomic center, 
      //  otherwise there will be many atomic center for the same type of AOs. eg. Zn(4s) in Mn:(ZnO)n.
      cout << setw(5) << AOCategory[i].center+1 << "-";
    else 
      cout << "      ";    //setw(6)
    
    cout << str_Atom[AOCategory[i].atomic] << "-" << str_L[AOCategory[i].l];
  };
  
  cout << endl;
  
  //print
  for(i=iStart;i<=iEnd;i++)
  {
    if( i == iHOMO ) cout << "---HOMO---LUMO----" << endl;   //print '---' below HOME.
    cout << setw(5) << i+1 ;
    cout << setw(8) << fixed << eigenEV[i];
    
    for(j=0;j<nAOTypes;j++) cout << setw(10) << fixed << CTperMO[i][j];
    
    cout << endl;
  }
  
};

//
void CCalcDOS::PrintCTperMO( const CCategoryAO* AOCategory )
{
  int i,j;
  cout.precision(5);
  
  //print head
  cout << "---iMO--Normalize--";
  if( nAOTypes < 12 )
  {
    for(j=0;j<nAOTypes;j++)  cout << setw(6) << AOCategory[j].atomic << "-" << AOCategory[j].l << " ";
  } else
  {
    for(j=0;j<6;j++)
      cout << setw(6) << AOCategory[j].atomic << "-" << AOCategory[j].l << " ";
    
    cout << " ...... ";
    
    for(j=nAOTypes-6;j<nAOTypes;j++)
      cout << setw(6) << AOCategory[j].atomic << "-" << AOCategory[j].l << " ";
    
  };
  cout << endl;
  
  //print CTperMO[][]
  for(i=0;i<nBasis;i++)
  {
    cout << setw(5)  << i+1 << " ";
    cout << setw(10) << fixed << normalCT[i] + 0.00001 << "     ";
    if( nAOTypes < 12 )
    {
      for(j=0;j<nAOTypes;j++) cout << setw(8) << fixed << CTperMO[i][j] << " ";
    } else
    {
      for(j=0;j<6;j++)  
        cout << setw(8) << fixed << CTperMO[i][j] << " ";
      
      cout << " ...... ";
      
      for(j=nAOTypes-6;j<nAOTypes;j++)
        cout << setw(8) << fixed << CTperMO[i][j] << " ";
    };
    
    cout << endl;
  };
  cout << "----------------------" << endl;
};

//
void CCalcDOS::FillEigen(const double *eigen)
{  
  for(int i=0; i<nBasis; i++) eigenEV[i] = au2ev * eigen[i]; 
};

//return eigenBins.size()
int CCalcDOS::FillEBins(double binStart, double binEnd, double binSize)
{
  eigenBins.clear();
  
  eigenNDOS.clear();
  eigenPDOS.clear();
  eigenCT.clear();
  
  eigenMOStart.clear();
  eigenMOEnd.clear();
  
  double binValue;
  int i=0;
  while(  ( binValue = ( binStart + binSize*(i++) )  ) < binEnd ) eigenBins.push_back( binValue );
  
  return i;
  
};

// DOS Total
void CCalcDOS::SumDOSTotal( const double* eigen, double binStart, double binEnd, double binSize )
{
  //init
  int i,j,k;
  FillEigen( eigen );
  int totalSteps = FillEBins( binStart, binEnd, binSize );
  
  int nState;
  int iAO=0;
  //find the first eigenEV[iAO] larger than eigenBins[0]
  while( eigenBins[0] > eigenEV[iAO] && iAO < nBasis ) iAO++; 
  
  for(i=0;i<totalSteps;i++)
  {
    //initialize
    nState = 0;
    eigenMOStart.push_back(iAO);
    
    //calc
    while( eigenBins[i] > eigenEV[iAO] && iAO < nBasis )
    {
      nState++;
      iAO++;
    };
    
    //assign
    eigenNDOS.push_back( nState );
    eigenMOEnd.push_back(iAO);
  };
};

//
void CCalcDOS::PrintDOSTotal(const CCategoryAO* AOCategory)
{
  cout << setw(8) << "E" << setw(8) << "NDOS" << endl;
  cout << "------------------------------------------" << endl;
  for(int i=1;i<eigenBins.size();i++)
    cout << setw(8) << setprecision(3) << fixed << eigenBins[i] << setw(8) << eigenNDOS[i] << endl;
  cout << "------------------------------------------" << endl;
};

//KB: start
//adding a separate routine for beta DOS, eigenvalues multiplied by -1.
void CCalcDOS::PrintDOSTotalBeta(const CCategoryAO* AOCategory)
{
  cout << setw(8) << "E" << setw(8) << "NDOS" << endl;
  cout << "------------------------------------------" << endl;
  for(int i=1;i<eigenBins.size();i++)
    cout << setw(8) << setprecision(3) << fixed << eigenBins[i] << setw(8) << -1 * eigenNDOS[i] << endl;
  cout << "------------------------------------------" << endl;
};
//KB: end

//DOS Partial
void CCalcDOS::SumDOSPartial( const double* eigen, double binStart, double binEnd, double binSize )
{
  //init
  int i,j,k;
  FillEigen( eigen );
  int totalSteps = FillEBins( binStart, binEnd, binSize );
  
  //
  int nState;
  double pState;
  vector<double>partial;
  
  int iAO=0;
  //find the first eigenEV[iAO] larger than eigenBins[0]
  while( eigenBins[0] > eigenEV[iAO] && iAO < nBasis ) iAO++; 
  
  for(i=0;i<totalSteps;i++)
  {
    //initialize
    partial.clear();
    nState = 0;
    pState = 0;
    for(j=0;j<nAOTypes;j++) partial.push_back( 0.0 );
    eigenMOStart.push_back( iAO );
    
    //calc
    while( eigenBins[i] > eigenEV[iAO] && iAO < nBasis )
    {
      nState++;
      for(j=0;j<nAOTypes;j++) partial[j] += CTperMO[iAO][j];
      iAO++; 
    };
    
    //assign
    eigenCT.push_back( partial );
    eigenNDOS.push_back( nState );
    eigenMOEnd.push_back( iAO );
    
    for(j=0;j<nAOTypes;j++) pState += partial[j];
    eigenPDOS.push_back( pState );
    
  };
  
};

//
void CCalcDOS::PrintDOSPartial( const CCategoryAO* AOCategory )
{
  int k,i,j;
  int format = 12 ;
  int totalBlocks = int ( double(nAOTypes) / double(format) + 0.99999 );
  
  for(k=0;k<totalBlocks;k++)
  {
    if( k > 0 ) cout << " *** Continuted **** "<< endl; 
    
    //print head
    //jwm cout << setw(8) << "E" << setw(6) << "Start" << setw(6) << "End" << setw(6) << "PDOS" << setw(6) << "NDOS";
	 cout << setw(8) << "Energy-a" << setw(8) << setw(8) << "total-a";
    
    for(j=k*format; j<(k+1)*format && j<nAOTypes; j++)  
       cout << setw(6) << str_Atom[AOCategory[j].atomic] << "-" << str_L[AOCategory[j].l] << "-a"; 
    
    cout << endl;
    
    //jwm cout << "-----------------------------------------------------------------------------------------------------------" << endl;
    
    //print eigenCT[][]
    for(i=1;i<eigenBins.size();i++)
    {
      cout << setw(8) << setprecision(3) << fixed << eigenBins[i] ;
      //jwm cout << setw(6) << eigenMOStart[i]+1;
      //jwm cout << setw(6) << eigenMOEnd[i]  +1;
      //jwm cout << setw(6) << setprecision(3) << fixed << eigenPDOS[i] ;
      cout << setw(6) << eigenNDOS[i] ;
      
      for(j=k*format; j<(k+1)*format && j<nAOTypes; j++) cout << setw(8) << setprecision(3) << fixed << eigenCT[i][j] ;
      
      cout << endl;
    };
  };
  //jwm cout << "-----------------------------------------------------------------------------------------------------------" << endl;
};


//KB start
//A copy of PrintDOSPartial for beta eigenvalues, multiplied by -1
void CCalcDOS::PrintDOSPartialBeta( const CCategoryAO* AOCategory )
{
  int k,i,j;
  int format = 12 ;
  int totalBlocks = int ( double(nAOTypes) / double(format) + 0.99999 );

  for(k=0;k<totalBlocks;k++)
  {
    if( k > 0 ) cout << " *** Continuted **** "<< endl;

    //print head
    //jwm cout << setw(8) << "E" << setw(8) << "Start" << setw(8) << "End" << setw(8) << "PDOS" << setw(8) << "NDOS";
	 cout << setw(8) << "Energy-b" << setw(8) << setw(8) << "total-b";

    for(j=k*format; j<(k+1)*format && j<nAOTypes; j++)
       cout << setw(8) << str_Atom[AOCategory[j].atomic] << "-" << str_L[AOCategory[j].l] << "-b";

    cout << endl;

    //jwm cout << "-----------------------------------------------------------------------------------------------------------" << endl;

    //print eigenCT[][]
    for(i=1;i<eigenBins.size();i++)
    {
      cout << setw(8) << setprecision(3) << fixed << eigenBins[i] ;
      //jwm cout << setw(8) << eigenMOStart[i]+1;
      //jwm cout << setw(8) << eigenMOEnd[i]  +1;
      //jwm cout << setw(8) << setprecision(3) << fixed << -1 * eigenPDOS[i] ;
      cout << setw(8) << -1 * eigenNDOS[i] ;

      for(j=k*format; j<(k+1)*format && j<nAOTypes; j++) cout << setw(8) << setprecision(3) << fixed << -1 * eigenCT[i][j] ;

      cout << endl;
    };
  };
  //jwm cout << "-----------------------------------------------------------------------------------------------------------" << endl;
};
//KB end

//
void CCalcDOS::SumElecDensity( int nElec )
{
  SumElecDensity( 0, nElec-1 );
};

//
void CCalcDOS::SumElecDensity( int iStart, int iEnd )
{
  int i,j;
  elecDensity.resize( nAOTypes );   //init
  for(i=0;i<nAOTypes;i++)
  {
    //sum CTperMO[][iAOTypes]
    for(j=iStart;j<=iEnd;j++)
    {
      elecDensity[i] += CTperMO[j][i];
    };
  };
};


//
void CCalcDOS::PrintElecDensity( const CCategoryAO* AOCategory )
{
  int i;
  //print head
  cout << "Total Electron Denisty::" << endl;
  for(i=0;i<nAOTypes;i++) 
  {
    if( AOCategory[i].center != -1 ) cout << setw(5) << AOCategory[i].center+1 << "-";   // only one atom center
    else cout << "      ";
    cout << str_Atom[AOCategory[i].atomic] << "-" << str_L[AOCategory[i].l];
  };
  cout << endl;
  
  //print elecDensity[];
  cout.precision(3);
  for(i=0;i<elecDensity.size();i++)  cout << setw(10) << fixed << elecDensity[i];
  cout << endl;
};
