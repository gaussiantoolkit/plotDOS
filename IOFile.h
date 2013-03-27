#ifndef INCLUDE_IOFILE
#define INCLUDE_IOFILE

#include "global.h"

class IOFile
{
public:
  bool OpenFile(string);                     //open a txt file;
  bool CloseFile();                          //close a file;
  
  void PrintFileInfo();                      //print out file information;
  
  bool SearchStr(string);                    //search the str;  
  int  PrintCurrentLine();                   //print out context till endl;
  
  bool ShiftCharsForward(unsigned int);      //forward txtfile read pointer;
  bool ShiftWordsForward(unsigned int);      //forward txtfile read pointer;
  bool ShiftLinesForward(unsigned int);      //forward txtfile read pointer;
  
protected:
  fstream inputFile;
  fstream outputFile;   //ToDo
  
  string  inputFileName;
  string outputFileName;
  int totalChars;
  int totalWords;       //ToDo
  int totalLines;
};

#endif
