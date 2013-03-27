#include "IOFile.h"

//open file
bool IOFile::OpenFile( string file )
{
  inputFileName = file;
  
  if( !inputFile.is_open() )
  {
    inputFile.open( inputFileName.c_str());
    if(inputFile.fail())
    {
      cerr << inputFileName << " doesn't exsit !." << endl;
      return false;
    };
  } 
  else
  {
    cerr << inputFileName << " has been opened already." << endl;
    return false;
  };  
  
  string tstr;
  int iLines = 0;
  
  while( !inputFile.eof()  ) 
  {
    getline(inputFile,tstr);            //getline();
    iLines++;
  };
  
  inputFile.clear();                    //reset
  totalChars = inputFile.tellg();       //count the '\n' in each line.
  totalLines = iLines-1; 
  inputFile.seekg(0);                   //reset
  
  return true;
};

//close file;
bool IOFile::CloseFile()
{
  inputFile.close();
  inputFileName = "";
  return true;
};

//print file info
void IOFile::PrintFileInfo()
{
  if( inputFile.is_open() )
  {
    cout << "'"<< inputFileName << "'" << " is opened!" << endl;
    cout << "Total Lines: " << setw(12) << totalLines << endl;
    cout << "Total Chars: " << setw(12) << totalChars << endl;
  } 
  else
  {
    cout << "There's no file opened!" << endl;
  };
};

//search string, there is no white-space in the matching string.
bool IOFile::SearchStr(string matching)
{
  string tstr ="";
  while( !inputFile.eof() && tstr != matching ) inputFile >> tstr;
  
  if( !inputFile.eof()) return true;
  else
  {
    //cerr <<"It's EOF. SearchStr() failed." << endl;
    inputFile.clear();
    return false;
  };
};

//print current line;
int IOFile::PrintCurrentLine()
{
  string tstr;
  int iCurrent = inputFile.tellg();
  
  getline(inputFile,tstr);
  cout << "{" << tstr << "}" << endl;
  
  if( !inputFile.eof() ) 
  {
    inputFile.seekg(iCurrent);     //reset
    return iCurrent;
  }    
  else
  {
    inputFile.clear();
    inputFile.seekg(iCurrent);    //reset
    return iCurrent;
  };
};

//shift forward by chars
bool IOFile::ShiftCharsForward(unsigned int n)   
{
  inputFile.seekg( int( inputFile.tellg() ) + n );
  
  if( !inputFile.eof() ) return true;
  else
  {
    cerr <<"It's EOF. ShiftCharsForward() failed." << endl;
    inputFile.clear();
    return false;
  };
};

//shift forward by words
bool IOFile::ShiftWordsForward(unsigned int n)   
{
  string tstr;
  while( ( n-- > 0 ) && ( !inputFile.eof() ) ) inputFile >> tstr;
  
  if( !inputFile.eof() ) return true;
  else
  {
    cerr <<"It's EOF. ShiftWordsForward() failed." << endl;
    inputFile.clear(); 
    return false;
  };
  
};

//shift forward by lines
bool IOFile::ShiftLinesForward(unsigned int n)
{
  string tstr;
  while( ( n-- > 0 ) && ( !inputFile.eof() ) ) getline(inputFile,tstr);
  
  if( !inputFile.eof()) return true;
  else 
  {
    cerr << "It's EOF. ShiftLinesForward() failed." << endl;
    inputFile.clear();
    return false;
  };
  
};



