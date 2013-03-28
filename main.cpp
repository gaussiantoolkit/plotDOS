#include <cstdlib>
#include "DOS.h"
#include "LOGFile.h"

// There is only one copy for 'MO-DOS' and 'MO-DOS-Diag'
// Please Comment 'DOS.cpp' Line:129, and UnComment Line:132 for no overlap matrix manipulation version

void DisplayHelp()
{
    cout<<" ***********************************************************************************"<<endl;
    cout<<" * Copyright(c) 2011 LiGroup@UWChem, All rights reserved.                          *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" * Last Updated by JWM on 2011-11-22                                               *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" * Name:                                                                           *"<<endl;
    cout<<" *   $plotDOS logfile.D2E total/partial alpha/beta (EVStart EVEnd StepSize)        *"<<endl;
    cout<<" *   $plotDOS logfile.D2E MO (nMO)                                                 *"<<endl;
    cout<<" *   $plotDOS logfile.D2E TM (nMO)                                                 *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" * Usage:                                                                          *"<<endl;
    cout<<" *    1.logfile must be created with 'GFInput IOP(3/33=1) POP(Full,NPA)' options   *"<<endl;
    cout<<" *      to include AO, overlap matrix and MO coefficent informations.              *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" *    2.In order to reduce the size of regular logfile. POP can be performed       *"<<endl;
    cout<<" *      with guess(read,only) in a seperated gauss job.                            *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" *    3.$sed 's/D+/E+/g;s/D-/E-/g;s/D 0/E00/' logfile > logfile.D2E                *"<<endl;
    cout<<" *      change all occurrences of 'D+/-' to 'E+/-'                                 *"<<endl;
    cout<<" *      change the first occurrence of 'D 0' to 'E00'                              *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" *    4.Limitation: N_ATOM_MAX = 1000; N_AO_MAX = 50000;                           *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" *    5.Don't include the following keywords in your com file.                     *"<<endl;
    cout<<" *            'NAtoms=' 'orientation:'  'Raffenetti'                               *"<<endl;
    cout<<" *            'basis' 'set' 'Overlap' 'Coefficients:'                              *"<<endl;
    cout<<" *                                                                                 *"<<endl;
    cout<<" *    6.If logfile readin geometry from checkpoint file, geom(check), then the     *"<<endl;
    cout<<" *      following line must be added to the logfile after the nuclear coordinates  *"<<endl;
    cout<<" *            'NAtoms= X'                                                          *"<<endl;
    cout<<" *      where 'X' is the number of atoms in your system.                           *"<<endl;
    cout<<" ***********************************************************************************"<<endl;
    exit(-1);
};

//** DOS() **//
void DOS( int ac, char** av, CLOGFileReadIn& logfileReadIn, CCataAOs& cataAOs)
{
  string strAlphaOrBeta, strJobType;
  double eigenStart, eigenEnd, stepsize;

  //parse DOS cmd line
  if( ac == 4 )
  {
    strJobType = *(av+2);
    strAlphaOrBeta = *(av+3);
    eigenStart = -15.0;
    eigenEnd   =  15.0;
    stepsize   =   0.2;
  } else if (ac == 7)
  {
    strJobType = *(av+2);
    strAlphaOrBeta = *(av+3);
    eigenStart = atof( *(av+4) );
    eigenEnd   = atof( *(av+5) );
    stepsize   = atof( *(av+6) );
    
    //check
    if( eigenStart > eigenEnd || stepsize < 0.001)
    {
      cout << "eigenStart, eigenEnd, stepsize ? " << endl;
      DisplayHelp();
      exit(-1);
    };
    
  } else 
  {
    DisplayHelp();
    exit(-1);
  };
  
  
  //Density of State
  CCalcDOS calcDOS( logfileReadIn.ReturnNBasis(), cataAOs.ReturnNAOTypes() );
  if( strAlphaOrBeta == "alpha" )
  {  
    if ( logfileReadIn.ReadInMOCoeffAlph() ) 
    {
      //logfileReadIn.PrintMOCoeffAlph();
      //logfileReadIn.PrintMOEigenAlph();
      ;
    } else
    {
      cerr << "Can't ReadIn MO Alpha Coeff !" << endl;
      DisplayHelp();
      exit(-1);
    };
    
    //Calc  
    if( strJobType == "total" )
    {
      calcDOS.SumDOSTotal( logfileReadIn.ReturnMOEigenAlph(), eigenStart, eigenEnd, stepsize );
      calcDOS.PrintDOSTotal( cataAOs.ReturnAOTypes() );
    } else // "partial"
    {
      calcDOS.CalcAOperMO( logfileReadIn.ReturnOverlap(), logfileReadIn.ReturnMOCoeffAlph() ) ;
      calcDOS.CalcNormalAO();
      //calcDOS.PrintAOperMO();
      
      calcDOS.SumCTperMO( cataAOs.ReturnAOTypes() );
      calcDOS.CalcNormalCT();
      //calcDOS.PrintCTperMO( cataAOs.ReturnAOTypes() );
      
      calcDOS.SumDOSPartial( logfileReadIn.ReturnMOEigenAlph(), eigenStart, eigenEnd, stepsize );
      calcDOS.PrintDOSPartial( cataAOs.ReturnAOTypes() );
    }
  }
  
  else if( strAlphaOrBeta == "beta" )
  
  {
    if ( logfileReadIn.ReadInMOCoeffBeta() ) 
    {
      //logfileReadIn.PrintMOCoeffBeta();
      //logfileReadIn.PrintMOEigenBeta();
      ;
    } else
    {
      cerr << "Can't ReadIn MO Beta Coeff !" << endl;
      DisplayHelp();
      exit(-1);
    };
    
    //Calc
    if( strJobType == "total" )
    {
      calcDOS.SumDOSTotal( logfileReadIn.ReturnMOEigenBeta(), eigenStart, eigenEnd, stepsize );
      calcDOS.PrintDOSTotalBeta( cataAOs.ReturnAOTypes() );
    } else // "partial"
    {
      calcDOS.CalcAOperMO( logfileReadIn.ReturnOverlap(), logfileReadIn.ReturnMOCoeffBeta() ) ;
      calcDOS.CalcNormalAO();
      //calcDOSAlph.PrintAOperMO();
      
      calcDOS.SumCTperMO( cataAOs.ReturnAOTypes() );
      calcDOS.CalcNormalCT();
      //calcDOSAlph.PrintCTperMO( cataAOs.ReturnAOTypes() );
      
      calcDOS.SumDOSPartial( logfileReadIn.ReturnMOEigenBeta(), eigenStart, eigenEnd, stepsize );
      calcDOS.PrintDOSPartialBeta( cataAOs.ReturnAOTypes() );
    } 
    
  } 
  
  else
  
  {
    cerr << "It's not a valid option. 'alpha or beta'? Please try again." << endl;
    DisplayHelp();
    exit(-1);
  };

};


//** ElecDensity() **//
void ElecDensity( int ac, char** av, CLOGFileReadIn& logfileReadIn, CCataAOs& cataAOs)
{
  //parze cmd line
  string strAlphaOrBeta;
  if( ac == 4 )
  {
    strAlphaOrBeta = *(av+3);
  } else 
  {
    DisplayHelp();
    exit(-1);
  };
  
  //total elec density
  CCalcDOS calcDOS( logfileReadIn.ReturnNBasis(), cataAOs.ReturnNAOTypes() );
  
  if( strAlphaOrBeta == "alpha" )
  {  
    if ( logfileReadIn.ReadInMOCoeffAlph() ) 
    {
      //Calc  
      calcDOS.CalcAOperMO( logfileReadIn.ReturnOverlap(), logfileReadIn.ReturnMOCoeffAlph() ) ;
      calcDOS.SumCTperMO( cataAOs.ReturnAOTypes() );
      calcDOS.SumElecDensity( logfileReadIn.ReturnNElecA() );
      calcDOS.PrintElecDensity( cataAOs.ReturnAOTypes() );
      
      //test
      //calcDOS.PrintAOperMO();
      
    } else
    {
      cerr << "Can't ReadIn MO Alpha Coeff !" << endl;
      DisplayHelp();
      exit(-1);
    };
    
  }
  
  else if( strAlphaOrBeta == "beta" )
  
  {
    if ( logfileReadIn.ReadInMOCoeffBeta() ) 
    {   
      //Calc
      calcDOS.CalcAOperMO( logfileReadIn.ReturnOverlap(), logfileReadIn.ReturnMOCoeffBeta() ) ;
      calcDOS.SumCTperMO( cataAOs.ReturnAOTypes() );
      calcDOS.SumElecDensity( logfileReadIn.ReturnNElecB() );
      calcDOS.PrintElecDensity( cataAOs.ReturnAOTypes() );
      
    } else
    {
      cerr << "Can't ReadIn MO Beta Coeff !" << endl;
      DisplayHelp();
      exit(-1);
    };
    
  } 
  
  else
  
  {
    cerr << "It's not a valid option. 'alpha or beta'? Please try again." << endl;
    DisplayHelp();
    exit(-1);
  };
  
};

//** TMCoeff()  **//
void TMCoeff( int ac, char** av, CLOGFileReadIn& logfileReadIn, CCataAOs& cataAOs )
{
  int nMO;
  int i, iMOStart, iMOEnd;
  
  //parse cmd line
  if( ac == 4 )
  {
    nMO = atoi( *(av+3) ); 
  } else
  {
    nMO = 1;   //default, interest MO is HOMO-1,HOMO & LUMO.
  };
  
  CCalcDOS calcDOS( logfileReadIn.ReturnNBasis(), cataAOs.ReturnNAOTypes() );
  
  
  if ( logfileReadIn.ReadInMOCoeffAlph() )       //alpha
  {
    cout << "***************alpha MOs*******************" << endl;
    if( nMO < ( logfileReadIn.ReturnNBasis() - logfileReadIn.ReturnNElecA() ) && nMO < logfileReadIn.ReturnNElecA() )
    {
      
      iMOStart = logfileReadIn.ReturnNElecA() - nMO - 1;
      iMOEnd   = logfileReadIn.ReturnNElecA() + nMO - 1;  
      
      for( i= iMOStart; i<=iMOEnd; i++)
      {
        calcDOS.CalcAOperMO( logfileReadIn.ReturnOverlap(), logfileReadIn.ReturnMOCoeffAlph(), i ) ;
        calcDOS.SumCTperMO( cataAOs.ReturnAOTypes(), i );
      };
      
      calcDOS.FillEigen( logfileReadIn.ReturnMOEigenAlph() );
      calcDOS.PrintCTperMO( cataAOs.ReturnAOTypes(), iMOStart, iMOEnd, logfileReadIn.ReturnNElecA() );
      
    } else
    {
      cerr << "It's not a valid nMO "<< endl;
      DisplayHelp();
      exit(-1);
    };
  };
  
  if ( logfileReadIn.ReadInMOCoeffBeta() )       //beta
  {
    cout << "***************beta MOs*******************" << endl;
    if( nMO < ( logfileReadIn.ReturnNBasis() - logfileReadIn.ReturnNElecB() ) && nMO < logfileReadIn.ReturnNElecB() )
    {
      
      iMOStart = logfileReadIn.ReturnNElecB() - nMO - 1;
      iMOEnd   = logfileReadIn.ReturnNElecB() + nMO - 1;  
      
      for( i= iMOStart; i<=iMOEnd; i++)
      {
        calcDOS.CalcAOperMO( logfileReadIn.ReturnOverlap(), logfileReadIn.ReturnMOCoeffBeta(), i ) ;
        calcDOS.SumCTperMO( cataAOs.ReturnAOTypes(), i );
      };
      
      calcDOS.FillEigen( logfileReadIn.ReturnMOEigenBeta() );
      calcDOS.PrintCTperMO( cataAOs.ReturnAOTypes(), iMOStart, iMOEnd, logfileReadIn.ReturnNElecB() );
      
    } else
    {
      cerr << "It's not a valid nMO "<< endl;
      DisplayHelp();
      exit(-1);
    };
  };
  

};

//** main() **//
int main( int argc, char** argv)
{
  if( argc < 3 ) 
  {
    DisplayHelp();
    exit(-1);
  };
  
  CLOGFileReadIn logfileReadIn;            //ReadIn
  
  if( logfileReadIn.OpenFile(*(argv+1)) 
      && logfileReadIn.ReadInAtoms()
      && logfileReadIn.ReadInBasisAO()
      && logfileReadIn.ReadInOverlap()  )
    {
      //logfileReadIn.PrintFileInfo();
      //logfileReadIn.PrintAtoms();
      //logfileReadIn.PrintBasisAO();  
      //logfileReadIn.PrintOverlap();
      ;
    } else 
    {
      DisplayHelp();
      exit(-1);
    };
  
  //Job Type
  string strJobType = *(argv+2);
  
  //Catagorize AOs
  CCataAOs cataAOs( logfileReadIn.ReturnNBasis() );
  
  if( strJobType == "total" || strJobType == "partial" ) 
  {
    cataAOs.CategorizeAO( logfileReadIn.ReturnAOBasis() );
    //cataAOs.PrintAOCategory();
    DOS( argc, argv, logfileReadIn,cataAOs);
  }
  else if ( strJobType == "density" )
  {
    cout << "nElecAlpha = " << logfileReadIn.ReturnNElecA() <<  endl;
    cout << "nElecBeta  = " << logfileReadIn.ReturnNElecB() <<  endl;
    
    cataAOs.CategorizeAO( logfileReadIn.ReturnAOBasis() );
    //cataAOs.PrintAOCategory();
    ElecDensity(argc, argv, logfileReadIn,cataAOs);
  }
  else if ( strJobType == "TM" )
  {
    cout << "nElecAlpha = " << logfileReadIn.ReturnNElecA() <<  endl;
    cout << "nElecBeta  = " << logfileReadIn.ReturnNElecB() <<  endl;
    
    cataAOs.CategorizeTM( logfileReadIn.ReturnAOBasis() );
    //cataAOs.PrintAOCategory();
    TMCoeff(argc, argv, logfileReadIn,cataAOs);
  }
  else if ( strJobType == "MO" ) 
  {
    cout << "nElecAlpha = " << logfileReadIn.ReturnNElecA() <<  endl;
    cout << "nElecBeta  = " << logfileReadIn.ReturnNElecB() <<  endl;
    cataAOs.CategorizeAO( logfileReadIn.ReturnAOBasis() );
    //cataAOs.PrintAOCategory();
    TMCoeff(argc, argv, logfileReadIn,cataAOs);
  }
  else 
  {
     cerr << "It's not a valid option. 'total, partial, TM, MO, density ... ' Please try again." << endl;
     DisplayHelp();
     exit(-1);
  };
  
};


