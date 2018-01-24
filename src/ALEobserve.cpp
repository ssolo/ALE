#include "ALE.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  cout << "ALEobserve using ALE v"<< ALE_VERSION <<endl;

  if (argc<2) 
    {
      cout << "usage:\n ./ALEobserve gene_tree_sample.newicks [burnin=0]" << endl;
      return 1;
    }

  string first_file=argv[1];
  boost::trim(first_file);	    
  string head=first_file;
  string ale_name=head+".ale";
  approx_posterior * ale;
  int burnin=0;  
  vector<string> ale_files;
  ale_files.push_back(first_file);
  for (int i=2;i<argc;i++)
    {
      string next_field=argv[i];
      vector <string> tokens;
      boost::split(tokens,next_field,boost::is_any_of("="),boost::token_compress_on);
      if (tokens[0]=="burnin")
	burnin=atoi(tokens[1].c_str());
      else
	ale_files.push_back(argv[i]);
    }
  ale=observe_ALE_from_file(ale_files,burnin);
  cout << "# observe "<< ale->observations << " tree(s) from: " <<  argv[1] ;
  for (int i=2;i<argc-1;i++)
    cout << " " << argv[i];
  cout << endl;
  cout << burnin<<" burn in per file discarded."<<endl;
  ale->save_state(ale_name);
  cout << "# saved in "<< ale_name<<endl;
  return 0;
}
