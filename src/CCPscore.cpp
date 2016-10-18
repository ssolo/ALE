
#include "exODT.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  string ale_file=argv[1];
  approx_posterior * ale;

  ale=load_ALE_from_file(ale_file);
  ifstream tree_stream (argv[2]);
  string tree;
  getline (tree_stream,tree);
  cout << ale->p(tree) << endl;

  return 1;
}
