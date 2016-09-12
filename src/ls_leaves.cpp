#include "ALE.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{

  map<string,int> names;
  for (int i=1;i<argc;i++)
    {
      ifstream file_stream1 (argv[i]);
      string tree;
      getline (file_stream1,tree);
      tree_type * T=TreeTemplateTools::parenthesisToTree(tree,false);
      vector <Node *> leaves=T->getLeaves();
      for (vector <Node *>::iterator it=leaves.begin();it!=leaves.end();it++)
	{
	  string name=(*it)->getName();
	  vector <string> tokens;
	  boost::split(tokens,name,boost::is_any_of("_"),boost::token_compress_on);
	  //name=tokens[0];

	  names[name]++;
	}      
    }
  for (map <string,int>::iterator it=names.begin();it!=names.end();it++)
    {
      cout << (*it).first << " " << (*it).second << endl;
    }
}
