#include "ALE.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  ifstream file_stream (argv[1]);
  string name;
  map <string,map<string,int> > keep;
  while(! file_stream.eof())
    {
      getline (file_stream,name);
      if (name.find("_")!=name.npos )
	{
	  vector <string> tokens;
	  boost::trim(name);	    
	  boost::split(tokens,name,boost::is_any_of(" \t"),boost::token_compress_on);
	  string sp=tokens[0];
	  string small_group=tokens[1];
	  string large_group=tokens[2];
	  keep[small_group][sp]=1;
	  keep[large_group][sp]=1;
	}
    
    }
  for (int i=2;i<argc;i++)
    {
      for ( map <string,map<string,int> >::iterator kit= keep.begin() ; kit!=keep.end() ; kit++ )
	{
	  ifstream file_stream1 (argv[i]);
	  string name=argv[i];
	  vector <string> tokens;
	  boost::split(tokens,name,boost::is_any_of("."),boost::token_compress_on);
	  //string outname="wol_hosts_"+tokens[0]+"_"+(*kit).first+".tree";
	  string outname="wol_paras_"+tokens[0]+"s_"+(*kit).first+".trees";
	  ofstream fout( outname.c_str() );

	  string tree;      
	  while(! file_stream1.eof())
	    {
	      getline (file_stream1,tree);
	      if (tree.find(")")!=tree.npos )
		{
		  
		  tree_type * T=TreeTemplateTools::parenthesisToTree(tree,false,"ID");
		  vector <Node *> leaves=T->getLeaves();	      		
		  for (vector <Node *>::iterator it=leaves.begin();it!=leaves.end();it++)
		    {
		      string name=(*it)->getName();
		      //cout << name << endl;
		      if (not (*kit).second.count(name)==1)
			{
			  //cout << "drop" << " " <<name << endl;
			  TreeTemplateTools::dropLeaf(*T,name);
			}
		    }    
		  fout << TreeTemplateTools::treeToParenthesis(*T) << endl;
  
		}
	    }
	}
    }

}
