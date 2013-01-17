#include "ALE.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  string ale_file=argv[1];
  string ale_name=ale_file+".ale";
  approx_posterior * ale;
  int burnin=0;
  if (argc>2) 
    burnin=atoi(argv[2]);
  ale=obsorve_ALE_from_file(ale_file,burnin);
  cout << "# observe "<< ale->observations << "trees from: " <<  argv[1] << endl;
  ale->save_state(ale_name);
  cout << "# saved in"<< ale_name<<endl;
  return 1;
}
