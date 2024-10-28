// This program is a helper program to parse command lines
// see program example_lineparser.C for an example of its use.

#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>  //uses exit

class lineparser
{
private:
  std::string pname;
  std::vector< char > names;
  std::vector<std::string> mytypes;
  std::vector<std::string> mydefaults;
  std::vector<std::string> mymean;

  public:
  lineparser() {;}
  ~lineparser() {;}

  void AddCommand(char c,std::string mtype, std::string mdef,
		  std::string meaning="")
  {
    names.push_back(c);
    mytypes.push_back(mtype);
    mydefaults.push_back(mdef);
    mymean.push_back(meaning);
  }

  bool Validate(char c, int &k)
  {
    bool doo = false;
    for (unsigned int j=0;j<names.size();j++)
      {
	if (c == names[j])
	  {
	    doo = true;
	    k = j;
	    break;
	  }
      }
    return doo;
  }

  bool Interpret(int argc, char **argv,int debug = 0)
  {
    if (debug >0) std::cerr << "Interpreting...." << std::endl;
    pname = argv[0];
    bool doo= false;
    for (int i=1; i< argc; i++)
      {
	char cu = *(argv[i]+1);
	if (debug>0) std::cerr << "Interpreting...." << i << " " 
			       << cu << " " << 
		       argv[i] << std::endl;       	
	int k;
	doo = Validate(cu,k);
	if ( doo ) 
	  {
	    mydefaults[k] = argv[i]+2;
	  }
	else
	  {
	    error(cu);
	  }
      }
    return doo;
  }

  void Print()
  {
    std::cout << "\t Commands Defined " << std::endl;
    std::cout << "Name  Type Default" << std::endl;
    
    for (unsigned int i=0;i<names.size();i++)
      {
	std::cout << names[i] << " " << mytypes[i] << " " << mydefaults[i]
		  << " " << mymean[i] << std::endl;
      }
  }

  bool GetBool(char c)
  {
    int index;
    bool ok = Validate(c,index);
    bool result = false;
    if (ok) 
      {
	if (mydefaults[index] == "true") result = true;
      }
    return result;
  }

  int GetInt(char c)
  {
    int index;
    bool ok = Validate(c,index);
    int result=0;
    if (ok) result= atoi(mydefaults[index].c_str());
    return result;
  }


  double GetDouble(char c)
  {
    int index;
    bool ok = Validate(c,index);
    double result=0;
    if (ok) result= atof(mydefaults[index].c_str());
    return result;

  }


  std::string GetString(char c)
  {
    int index;
    bool ok = Validate(c,index);
    std::string result="";
    if (ok) result= mydefaults[index];
    return result;
  }

void error(char c)
{
  std::cerr << "Error. Command  " << c
	    << " not defined" << std::endl;

  std::cerr << "Usage: "<<std::endl;
  std::cerr << "\t " << pname << " [flags] "<<std::endl;
  std::cerr << " \t where flags can be: " << std::endl;
  
  for (unsigned int i=0;i< names.size(); i++)
    {
      std::cerr << "\t\t -" << names[i] << " " << mymean[i] << "  Defined Value: " << mydefaults[i]
		<< std::endl;
    }
  exit(1);
}

};
