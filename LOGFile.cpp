#include "LOGFile.h"

//** ReadIn Atoms Info **//
bool CLOGFileReadIn::ReadInAtoms()
{
  int i;
  int iCurrent = inputFile.tellg();
  inputFile.seekg(0);
  
  //nAtoms
  if( SearchStr("NAtoms=") )
  {
    //PrintCurrentLine();
    inputFile >> nAtoms;
    //PrintCurrentLine();
  } else
  {
    cerr << "Can't find 'NAtoms=' in 'logfile'! " << endl;
    return false;
  };
  
  //atom coordinates
  if( SearchStr("orientation:") )
  {
    //PrintCurrentLine();
    ShiftLinesForward(5);
    for(i=0;i<nAtoms;i++)
    { 
      ShiftWordsForward(1);
      atoms[i].center=i;    //start from '0'
      inputFile >> atoms[i].atomic >> atoms[i].atomic_type;
      inputFile >> atoms[i].xyz[0] >> atoms[i].xyz[1] >> atoms[i].xyz[2]; 
      ShiftLinesForward(1);
    };
   // PrintCurrentLine();
  } else
  {
    cerr << "Can't find 'orientation:' in 'logfile'! " << endl;
    return false;
  };
  
  //nAOBasis,nElecAlpha,nElecBeta
  if( SearchStr("Raffenetti") )
  {
    //PrintCurrentLine(); 
    ShiftLinesForward(2);  
    inputFile >> nAOBasis;
    
    ShiftLinesForward(1); 
    inputFile >> nElecA;
    nVirtA = nAOBasis - nElecA;
    
    ShiftWordsForward(2); 
    inputFile >> nElecB;
    nVirtB = nAOBasis - nElecB;
    //PrintCurrentLine();
  }else
  {
    cerr << "Can't find 'Raffenetti' in 'logfile'! " << endl;
    return false;
  };
  //reset the original file ptr.
  inputFile.seekg(iCurrent);
  return true;
};

//print
void CLOGFileReadIn::PrintAtoms()
{
  cout << "nAtoms:      " << nAtoms   << endl;
  cout << "nAOBasis:    " << nAOBasis << endl;
  cout << "nElecAlpha:  " << nElecA   << endl;
  cout << "nElecBeta:   " << nElecB   << endl;
  
  //xyz
  cout << "***Atoms***" << endl;
  cout << "Center\t" << "Atomic\t"<< "     X    \t" << "    Y     \t" << "    Z     " << endl;
  cout.precision(5);
  for(int i=0;i<nAtoms;i++)
  {
    cout << setw(4)  << atoms[i].center+1 << "\t" << setw(4) << atoms[i].atomic << "\t";
    cout << setw(10) << fixed << atoms[i].xyz[0] << "\t";
    cout << setw(10) << fixed << atoms[i].xyz[1] << "\t";
    cout << setw(10) << fixed << atoms[i].xyz[2] << endl;
  };
  
};

//** ReadIn Basis Set Info **//
bool CLOGFileReadIn::ReadInBasisAO()
{
  int iAO = 0;
  int i,j;
  int iCurrent = inputFile.tellg();
  inputFile.seekg(0);
  
  if( SearchStr("basis") && SearchStr("set") )
  {
    int nSTO;
    int iAtoms;
    string tSPDF="";
    //PrintCurrentLine();
    ShiftLinesForward(1);
    for(i=0;i<nAtoms;i++)
    { 
      inputFile >> iAtoms;    //iAtoms == i+1;
      if( iAtoms != i+1 )
      {
        cerr << "ReadIn AO Error! iAtoms:" << iAtoms << "\t"<< "i:" << i <<  endl; 
        return false;
      };
      
      ShiftLinesForward(1);
      
      while( tSPDF != "****" )
      {
        inputFile >> tSPDF;
        
        if( tSPDF != "****" ) 
        {              
          inputFile >> nSTO;
          
          //S P SP or D, ToDO 'F'
          if( tSPDF == "S")
          {
            for(j=0;j<1;j++)
            {
              AOBasis[iAO].center = i;
              AOBasis[iAO].atomic = atoms[i].atomic;
              AOBasis[iAO].l = 0;
              iAO++;
            };
          } 
          else if (tSPDF == "P" )
          {
            for(j=0;j<3;j++)
            {
              AOBasis[iAO].center = i;
              AOBasis[iAO].atomic = atoms[i].atomic;
              AOBasis[iAO].l = 1;
              iAO++;
            };
          } 
          else if (tSPDF == "SP")
          {
            for(j=0;j<1;j++)
            {
              AOBasis[iAO].center = i;
              AOBasis[iAO].atomic = atoms[i].atomic;
              AOBasis[iAO].l = 0;
              iAO++;
            };
            for(j=0;j<3;j++)
            {
              AOBasis[iAO].center = i;
              AOBasis[iAO].atomic = atoms[i].atomic;
              AOBasis[iAO].l = 1;
              iAO++;
            };            
          } 
          else if (tSPDF == "D" )
          {
            for(j=0;j<5;j++)
            {
              AOBasis[iAO].center = i;
              AOBasis[iAO].atomic = atoms[i].atomic;
              AOBasis[iAO].l = 2;
              iAO++;
            };
          } 
          else 
          {
            cout << "ReadIn AO Error!" << endl; 
            return false;
          };              
        }
        else nSTO = 0;            
        ShiftLinesForward(nSTO+1);
      };
      tSPDF = "";
    };
  }
  else
  {
    cerr << "Can't find 'basis' in 'logfile'! " << endl;
    return false;
  };
  
  //PrintCurrentLine();
  inputFile.seekg(iCurrent);
  return true;
};

//print
void CLOGFileReadIn::PrintBasisAO()
{
  cout << "***AO Basis***" << endl;
  cout << "iBasis\t" << "Center\t" <<"Atomic\t" << "L" << endl;
  for(int i=0; i<nAOBasis; i++)
  {
    cout << setw(4) << i+1 << "\t";
    cout << setw(4) << AOBasis[i].center+1 << "\t";
    cout << setw(4) << AOBasis[i].atomic   << "\t";
    cout << setw(1) << AOBasis[i].l << endl;
  };
};

//** ReadIn Overlap Matrix **//
bool CLOGFileReadIn::ReadInOverlap()
{
  //jwm cout << "ReadIn Overlap Matrix :: ";
  int iCurrent = inputFile.tellg();
  inputFile.seekg(0);
  
  string tstr;
  int i,j;
  int nFormat = 5;                // '5' is format of logfile
  int nBlock, nFrags, iBlock;
  nBlock = nAOBasis / nFormat ;
  nFrags = nAOBasis % nFormat ;

  int iDotStep = nBlock / 20 + 1;  //while doing time consuming jobs, '20' dots will be display one by one.
  
  if( SearchStr( "Overlap" ) )
  {
    //PrintCurrentLine();
    //readin overlap matrix step 1: Blocks
    for(iBlock=0;iBlock<nBlock;iBlock++)
    {
      if( iBlock % iDotStep == 0 ); //jwm cout << ". " << flush;
      ShiftLinesForward(2);      
      //first 'nFormat' lines
      for( i=iBlock*nFormat; i<(iBlock+1)*nFormat; i++)
      {
        ShiftWordsForward(1);
        for(j=iBlock*nFormat;j<=i;j++)                 
          inputFile >> overlapMatrix[i][j];                  
      };
      
      //the remaining lines
      for( ; i<nAOBasis; i++)
      {
        ShiftWordsForward(1); 
        for(j=iBlock*nFormat;j<(iBlock+1)*nFormat;j++)
          inputFile >> overlapMatrix[i][j];
      };            
    };
    
    //jwm cout << "." <<endl;
    //readin overlap matrix step 2: the Frags Block   
    if(nFrags>0){
      ShiftLinesForward(2);
      for( i = iBlock*nFormat;i<(iBlock*nFormat + nFrags);i++)
      {
        ShiftWordsForward(1);
        for(j=iBlock*nFormat;j<=i;j++)
          inputFile >> overlapMatrix[i][j];
      };
    };
    
    //fill overlapMatrix[i][j] when i<j
    for( i=0; i<nAOBasis; i++)
      for( j=0; j<i; j++)
      overlapMatrix[j][i] = overlapMatrix[i][j]; 
    
    ShiftLinesForward(1);
    
  }  else
  {
    cerr << "Can't find 'Overlap' in 'logfile'! " << endl;
    return false;
  };
  
  //PrintCurrentLine();
  inputFile.seekg(iCurrent);
  return true;
};

//print
void CLOGFileReadIn::PrintOverlap()
{
  cout << "*** Overlap Matrix, if nAOBasis > 12, only first and last six AO are printed! ***" << endl;
  int i,j;
  
  if( nAOBasis < 12 )
  {
    //half matrix
    for(i=0;i<nAOBasis;i++)
    {
      cout << setw(6) << i+1 << "\t";
      for(j=0;j<=i;j++)
        cout << setw(15) << scientific << overlapMatrix[i][j] << "\t";
      cout << endl;
    };
    
    cout << endl;
    
    //full matrix
    for(i=0;i<nAOBasis;i++)
    {
      cout << setw(6) << i+1 << "\t";
      for(j=0;j<nAOBasis;j++)
        cout << setw(15) << scientific << overlapMatrix[i][j] << "\t";
      cout << endl;
    };
  } else
  {
    //nAOBasis > 12
    for(i=0;i<nAOBasis;i++)
    {
      cout << setw(6) << i+1 << "\t";
      for(j=0;         j<6;       j++) cout << setw(15) << scientific << overlapMatrix[i][j] << "\t";
      cout << "  ......\t"; 
      for(j=nAOBasis-6;j<nAOBasis;j++) cout << setw(15) << scientific << overlapMatrix[i][j] << "\t";
      cout << endl;
    };
    
  };
};

//after Calling SearchStr( "Coefficients:" ), then ...
void CLOGFileReadIn::ReadInMOCoeff( double (*coeff)[nBASIS_MAX], double* eigen )
{
  int i,j;
  int nFormat = 5;           // '5' is format of logfile
  int nBlock, nFrags, iBlock;
  nBlock = nAOBasis / nFormat ;
  nFrags = nAOBasis % nFormat ;
  string tstr;
  
  int iDotStep = nBlock / 20 + 1;  //while doing time consuming jobs, '20' dots will be display one by one
  
  ShiftLinesForward(1);
  
  //readin coeff matrix step 1: Blocks
  for(iBlock=0;iBlock<nBlock;iBlock++)
  { 
    if( iBlock % iDotStep == 0 ); //jwm cout << ". " << flush;
    
    ShiftLinesForward(2);
    ShiftWordsForward(2);
    
    for(i=0;i<nFormat;i++) inputFile >> eigen[iBlock*nFormat+i];
    
    ShiftLinesForward(1);      
    for(j=0;j<nAOBasis;j++)
    {
      ShiftCharsForward(20);   //format of logfile
      for(i=0;i<nFormat;i++) inputFile >> coeff[iBlock*nFormat+i][j];
      ShiftLinesForward(1);
    };    	  	
  };
  
  //jwm cout << "." << endl;
  
  //readin coeff matrix step2: Frags Block
  if(nFrags>0)
  {
    ShiftLinesForward(2);
    ShiftWordsForward(2);
    
    for(i=0;i<nFrags;i++) inputFile >> eigen[iBlock*nFormat+i];
    
    ShiftLinesForward(1);     //end of line
    for(j=0;j<nAOBasis;j++)
    {
      ShiftCharsForward(20);  //format of logfile
      for(i=0;i<nFrags;i++) inputFile >> coeff[iBlock*nFormat+i][j];
      ShiftLinesForward(1);   //end of line
    }
  };
  
};

// ReadIn Alpha MO Coeff
bool CLOGFileReadIn::ReadInMOCoeffAlph()
{
  int iCurrent = inputFile.tellg();
  inputFile.seekg(0);
  
  if( SearchStr("Coefficients:") )
  {
    //PrintCurrentLine();
    //jwm cout << "ReadIn alpha MO coeff :: ";
    ReadInMOCoeff( coeffMOAlph, eigenAlph );
  } else
  {
    cerr << "Can't find 'Coefficients:' in 'logfile'! " << endl;
    return false;
  };
  
  //PrintCurrentLine();
  inputFile.seekg(iCurrent);
  return true;
};

// ReadIn Beta MO Coeff
bool CLOGFileReadIn::ReadInMOCoeffBeta()
{
  int iCurrent = inputFile.tellg();
  inputFile.seekg(0);
  
  if( SearchStr("Coefficients:") )
  {
    
    if( SearchStr("Coefficients:") ) 
    { 
      //PrintCurrentLine();
      //jwm cout << "ReadIn beta MO coeff  :: ";
      ReadInMOCoeff( coeffMOBeta, eigenBeta );
    }else
    {
      cerr << "Can't find 'Coefficients:' twice in 'logfile'! " << endl;
      return false;
    };
    
  } else
  {
    cerr << "Can't find 'Coefficients:' in 'logfile'! " << endl;
    return false;
  };
  
  //PrintCurrentLine();
  inputFile.seekg(iCurrent);
  return true;
};
  
//print
void CLOGFileReadIn::PrintMOCoeffAlph()
{
  int i,j;
  cout.precision(6);
  cout << "*** MO Coeff Alplha, if nAOBasis > 12, only first and last six MO are printed! ***" << endl ;
  for(i=0;i<nAOBasis;i++)
  {
    cout << setw(6) << i+1 << "\t";
    if( nAOBasis < 12 )
    {
      for(j=0;j<nAOBasis;j++)  cout << setw(10) << fixed << coeffMOAlph[i][j] << "\t";
    } else
    {
      for(j=0;         j<6;       j++)  cout << setw(10) << fixed << coeffMOAlph[i][j] << "\t";
      cout << "  ......\t";
      for(j=nAOBasis-6;j<nAOBasis;j++)  cout << setw(10) << fixed << coeffMOAlph[i][j] << "\t";
    };
    cout << endl;
  };
};

//print
void CLOGFileReadIn::PrintMOCoeffBeta()
{
  int i,j;
  cout.precision(6);
  cout << "*** MO Coeff Beta, if nAOBasis > 12, only first and last six MO are printed! ***" << endl ;
  for(i=0;i<nAOBasis;i++)
  {
    cout << setw(6) << i+1 << "\t";
    if( nAOBasis < 12 )
    {
      for(j=0;j<nAOBasis;j++)  cout << setw(10) << fixed << coeffMOBeta[i][j] << "\t";
    } else
    {
      for(j=0;         j<6;       j++)  cout << setw(10) << fixed << coeffMOBeta[i][j] << "\t";
      cout << "  ......\t";
      for(j=nAOBasis-6;j<nAOBasis;j++)  cout << setw(10) << fixed << coeffMOBeta[i][j] << "\t";
    };
    cout << endl;
  };
};

//print
void CLOGFileReadIn::PrintMOEigenAlph()
{
  int i;
  cout.precision(6);
  cout << "***  Eigen Value Alpha(au) ***" << endl ;
  for(i=0;i<nAOBasis;i++)
  {
    cout << setw(6) << i+1 << "\t" ;
    cout << setw(10) << fixed << eigenAlph[i] << endl;
  };
};

//print
void CLOGFileReadIn::PrintMOEigenBeta()
{
  int i;
  cout.precision(6);
  cout << "***  Eigen Value Beta(au) ***" << endl ;
  for(i=0;i<nAOBasis;i++)
  {
    cout << setw(6) << i+1 << "\t" ;
    cout << setw(10) << fixed << eigenBeta[i] << endl;
  };
};


