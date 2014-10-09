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
  vector <string> tokens;
  boost::trim(first_file);	    
  boost::split(tokens,first_file,boost::is_any_of("."),boost::token_compress_on);
  string head=tokens[0];
  
  vector<string> ale_files;
  ale_files.push_back(first_file);
  for (int i=2;i<argc-1;i++)
    ale_files.push_back(argv[i]);
  string ale_name=head+".ale";
  approx_posterior * ale;
  int burnin=500;
  if (argc>2) 
    burnin=atoi(argv[argc-1]);
  ale=observe_ALE_from_file(ale_files,burnin);
  cout << "# observe "<< ale->observations << "trees from: " <<  argv[1] ;
  for (int i=2;i<argc-1;i++)
    cout << " " << argv[i];
  cout << endl;
  ale->save_state(ale_name);
  cout << "# saved in "<< ale_name<<endl;
  return 1;
}
