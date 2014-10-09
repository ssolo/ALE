
#include "exODT.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  string ale_file=argv[1];
  approx_posterior * ale;

  ale=load_ALE_from_file(ale_file);
  //string G="((((((((A,B),C),D),E),F),G),J),H);";

  cout << ale->count_trees()<< endl;
  //cout << ale->count_all_trees(ale->Gamma)<< endl;

  return 1;
}
/*
{
  //string ale_file=argv[2];
  approx_posterior * ale;
  approx_posterior * ale2;

  //ale=observe_ALE_from_file(ale_file);
  string G="((((((A,B),C),D),E),F),G);";
  //G="((A,B),(C,D));";
  ale=observe_ALE_from_string(G);
  vector <string> all_trees=ale->all_trees(ale->Gamma);

  ale2=observe_ALE_from_strings(all_trees);

  cout << ale->count_trees()<< endl;
  cout << ale->count_all_trees(ale->Gamma)<< endl;

  cout << ale2->count_trees()<< endl;
  cout << ale2->count_all_trees(ale->Gamma)<< endl;

  return 1;
}
*/
