#include "ALE.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  cout << "ALEadd using ALE v"<< ALE_VERSION <<endl;

  if (argc<2) 
    {
      cout << "usage:\n ./ALEadd ale_file.ale gene_tree_sample.newicks [weight=1] [burnin=0] [every=1] [until=end] [outfile=filename]" << endl;
      return 1;
    }

  string ale_file=argv[1];
  string trees_file=argv[2];
  boost::trim(ale_file);	    
  boost::trim(trees_file);	    
  string ale_name=ale_file;
  approx_posterior * ale;
  int burnin=0,every=1,until=-1;
  scalar_type weight = 1;
  
  for (int i=3;i<argc;i++)
    {
      string next_field=argv[i];
      vector <string> tokens;
      boost::split(tokens,next_field,boost::is_any_of("="),boost::token_compress_on);
      if (tokens[0]=="burnin")
	burnin=atoi(tokens[1].c_str());
      else if (tokens[0]=="every")
	every=atoi(tokens[1].c_str());
      else if (tokens[0]=="until")
	until=atoi(tokens[1].c_str());
      else if (tokens[0]=="weight")
	weight=atof(tokens[1].c_str());
      else if (tokens[0]=="outfile")
	ale_name=tokens[1];
    }
  ale=load_ALE_from_file(ale_file);
  cout <<  "." << endl;


  vector<string> trees;
  ifstream file_stream (trees_file.c_str());
  int tree_i=0;  
  if (file_stream.is_open())  //  ########## read trees ############
    {
      while (! file_stream.eof())
	{
	  string line;
	  getline (file_stream,line);
	  if (line.find("(")!=line.npos)
	    {
	      if (tree_i>=burnin and tree_i%every==0)
		{
		  cout << line;
		  trees.push_back(line);		 
		}
		tree_i++;
	    }
	}
    }

  cout <<  ".." << endl;

  vector<string> observe_trees;
  if (until==-1)
    until=trees.size();
  for (int i=0;i<min((int)trees.size(),until);i++)
    {
      observe_trees.push_back(trees[i]);
    }
  ale->observation(observe_trees,false,weight);


  cout <<"# " << observe_trees.size() <<  " new tree(s) observed with weight "<<weight<<" from: " <<  argv[2] ;
  cout <<"; " << burnin<<" trees burnin discarded."<<endl;
  
  cout << "# .ale with "<< ale->observations << " tree(s) from: " <<  argv[1] << " and " << argv[2] << endl;
  ale->save_state(ale_name);
  cout << "# saved in "<< ale_name<<endl;
  
  trees.clear();
  observe_trees.clear();

  return 1;
}
