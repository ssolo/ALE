#include "ALE_util.h"
using namespace std;
using namespace bpp;

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  if (ifile) return true;
    return false;
};

approx_posterior * observe_ALE_from_file(vector<string> fnames, int burnin,int every,int until)
{


  vector<string> trees;
  for (vector<string>::iterator it=fnames.begin();it!=fnames.end();it++)
    {
      string fname=(*it);
      ifstream file_stream (fname.c_str());
      int tree_i=0;  
      if (file_stream.is_open())  //  ########## read trees ############
	{
	  while (! file_stream.eof())
	    {
	      string line;
	      getline (file_stream,line);
	      if (line.find("(")!=line.npos )
		{
		  tree_i++;
		  if (tree_i>burnin and tree_i%every==0) trees.push_back(line);			     
		}
	    }
	}
    }
  if (trees.size()<1)
    return NULL;

  vector<string> observe_trees;
  if (until==-1)
    until=trees.size();
  for (int i=0;i<min((int)trees.size(),until);i++)
    observe_trees.push_back(trees[i]);
  cout << observe_trees.size() << endl;

  approx_posterior* ale=new approx_posterior(trees[0]);// NO del-loc

  ale->observation(observe_trees);

  trees.clear();
  observe_trees.clear();
  
  return ale;     
}
 
approx_posterior * observe_ALE_from_file(string fname, int burnin,int every,int until)
{

  vector<string> trees;
  ifstream file_stream (fname.c_str());
  int tree_i=0;  
  if (file_stream.is_open())  //  ########## read trees ############
    {
      while (! file_stream.eof())
	{
	  string line;
	  getline (file_stream,line);
	  if (line.find("(")!=line.npos)
	    {
	      tree_i++;
	      if (tree_i>burnin and tree_i%every==0) trees.push_back(line);			     
	    }
	}
    }

  approx_posterior* ale=new approx_posterior(trees[0]);// NO del-loc

  vector<string> observe_trees;
  if (until==-1)
    until=trees.size();
  for (int i=0;i<min((int)trees.size(),until);i++)
    observe_trees.push_back(trees[i]);

  ale->observation(observe_trees);

  trees.clear();
  observe_trees.clear();

  return ale;     
}

approx_posterior * observe_ALE_from_string(string tree)
{
  vector<string> trees;
  trees.push_back(tree);
  approx_posterior* ale=new approx_posterior(trees[0]);// NO del-loc  
  ale->observation(trees,false);
  return ale;     
}
approx_posterior * observe_ALE_from_strings(vector<string> trees)
{
  approx_posterior* ale=new approx_posterior(trees[0]);// NO del-loc  
  ale->observation(trees,false);
  return ale;     
}

//this function is specific to a dataset with a particular problem with species name seporators 
//DO NOT USE AS IS
approx_posterior * observe_ALE_from_nexus(string fname, int burnin,int every,int until)
{
  vector<string> trees;
  map <string,string> translate;
  ifstream file_stream (fname.c_str());
  int tree_i=0;  
  cout << "reading nexus." <<endl;

  if (file_stream.is_open())  //  ########## read trees ############
    {
      bool header=true;
      while (! file_stream.eof())
	{
	  string line;
	  getline (file_stream,line);
	  if (line.find("tree gen")!=line.npos)
	    header=false;
	  if (header)
	    {
	      if (line.find("Param:")==line.npos and line.find("begin trees")==line.npos and line.find("translate")==line.npos and line.find("#NEXUS")==line.npos and line.find("ID")==line.npos)
		{
		  vector <string> tokens;
		  boost::trim(line);	    
		  boost::split(tokens,line,boost::is_any_of(",; "),boost::token_compress_on);
		  vector <string> name_tokens;
		  //cout << tokens[1] << endl;;
		  boost::split(name_tokens,tokens[1],boost::is_any_of("_"),boost::token_compress_on);
		  string new_name=name_tokens[0]+"-"+name_tokens[1];
		  for (int i = 2 ; i<(int)name_tokens.size();i++)
		      new_name+="_"+name_tokens[i];
		  translate[tokens[0]]=new_name;
		}
	    }
	  else
	    {
	      if (line.find("(")!=line.npos)
		{
		  vector <string> tokens;
		  boost::trim(line);	    
		  boost::split(tokens,line,boost::is_any_of(" "),boost::token_compress_on);
		  tree_type * tree=TreeTemplateTools::parenthesisToTree(tokens[4],false);
		  vector <Node*> leaves=tree->getLeaves();
		  for (vector <Node*> :: iterator it=leaves.begin();it!=leaves.end();it++)
		    (*it)->setName(translate[(*it)->getName()]);
		  tree_i++;
		  if (tree_i>burnin and tree_i%every==0) trees.push_back(TreeTemplateTools::treeToParenthesis(*tree));
		  delete tree;
		  
		}
	    }
	}
    }
  cout << "translated nexus." <<endl;
  approx_posterior* ale=new approx_posterior(trees[0]);// NO del-loc
  vector<string> observe_trees;
  if (until==-1)
    until=trees.size();
  for (int i=0;i<min((int)trees.size(),until);i++)
    observe_trees.push_back(trees[i]);

  cout << "start observe." <<endl;
  ale->observation(observe_trees,false);
  trees.clear();
  observe_trees.clear();
  return ale;     
}

approx_posterior * load_ALE_from_file(string fname)
{
  approx_posterior* ale=new approx_posterior();// NO del-loc
  ale->load_state(fname);
  return ale;     
}

string save_ALE_to_file(string fname)
{
  vector<string> trees;
  ifstream file_stream (fname.c_str());
  
  if (file_stream.is_open())  //  ########## read trees ############
    {
      while (! file_stream.eof())
	{
	  string line;
	  getline (file_stream,line);
	  if (line.find("(")!=line.npos)
	    trees.push_back(line);			     
	}
    }
  approx_posterior* ale=new approx_posterior(trees[0]);// del-loc
  ale->observation(trees,false);
  
  vector <string> tokens;
  boost::trim(fname);	    
  boost::split(tokens,fname,boost::is_any_of("."),boost::token_compress_on);
  fname=tokens[0];

  ofstream fout( (fname+".trees").c_str() );
  fout << "#tree" << " " << "pp" << " " << "alepp" << endl;
  
  for ( map< string , int >::iterator it=ale->tree_counts.begin();it!=ale->tree_counts.end();it++)
    fout << (*it).first << " " << (*it).second/ale->observations << " " << ale->p((*it).first) << endl;
  
  ale->save_state(fname+".ale");  
  delete ale;
  return fname+".ale";
}

string canonical_branch_lengths( string Sstring)
{
  tree_type * S =TreeTemplateTools::parenthesisToTree(Sstring,false);//del-loc
  vector <Node *> nodes=S->getNodes();//del-loc
  map <Node*,scalar_type> node2height;//del-loc
  map <scalar_type,Node*> height2node;//del-loc
  for (vector <Node *>::iterator it=nodes.begin();it!=nodes.end();it++)
    {      
      if ((*it)->isLeaf())
	{
	  node2height[(*it)]=0;	  
	}
      else
	{	  
	  vector<Node*> sons=(*it)->getSons();
	  scalar_type h0 = sons[0]->getDistanceToFather() + node2height[sons[0]];
	  scalar_type h1 = sons[1]->getDistanceToFather() + node2height[sons[1]];
	  if (abs(h0-h1)>1e-3)
	    {
	      cout << " tree is not ultrametric! with diff " << abs(h0-h1)<< endl;
	      h0=max(h0,h1);
	    }
	  node2height[(*it)]=h0;	  
	  if (height2node.count(h0))
	    {
	      cout << " tree is degenerate! at height " << h0 <<endl;
	      height2node[h0+1e-6]=(*it);	      
	    }
	  else
	    height2node[h0]=(*it);	      	  
	}
    }
  //map <Node*,int> node2rank;//del-loc
  //map <int,Node*> rank2node;//del-loc
  int rank=0;

  map <Node*,scalar_type> new_height;//del-loc
  int n=S->getNumberOfLeaves();
  scalar_type rank_height=0;
  
  for (map <scalar_type,Node *>::iterator hit=height2node.begin();hit!=height2node.end();hit++)
    {
      rank_height+=1.0/(scalar_type)(n-rank);
      rank+=1;
      //rank2node[rank]=(*hit).second;
      //node2rank[(*hit).second]=rank;
      new_height[(*hit).second]=rank_height;
    }
  
  for (vector <Node *>::iterator it=nodes.begin();it!=nodes.end();it++)
    if ((*it)->hasFather())
      {
      if ((*it)->isLeaf())
	{
	  scalar_type fathers_height= new_height[ (*it)->getFather() ];
	  (*it)->setDistanceToFather(fathers_height/rank_height);	  
	}
      else
	{
	  scalar_type fathers_height= new_height[ (*it)->getFather() ];
	  (*it)->setDistanceToFather((fathers_height- new_height[ (*it)])/rank_height);	  
	}
      }
  node2height.clear();
  height2node.clear();

  //node2rank.clear();
  //rank2node.clear();

  new_height.clear();

  nodes.clear();
  Sstring=TreeTemplateTools::treeToParenthesis(*S);
  delete S;

  return Sstring;
}



void canonical_branch_lengths( tree_type * S )
{
	vector <Node *> nodes=S->getNodes();//del-loc
	map <Node*,scalar_type> node2height;//del-loc
	map <scalar_type,Node*> height2node;//del-loc
	for (vector <Node *>::iterator it=nodes.begin();it!=nodes.end();it++)
    {
		if ((*it)->isLeaf())
		{
			node2height[(*it)]=0;
		}
		else
		{
			vector<Node*> sons=(*it)->getSons();
			scalar_type h0 = sons[0]->getDistanceToFather() + node2height[sons[0]];
			scalar_type h1 = sons[1]->getDistanceToFather() + node2height[sons[1]];
			if (abs(h0-h1)>1e-3)
			{
				cout << " tree is not ultrametric! with diff " << abs(h0-h1)<< endl;
				h0=max(h0,h1);
			}
			node2height[(*it)]=h0;
			if (height2node.count(h0))
			{
				cout << " tree is degenerate! at height " << h0 <<endl;
				height2node[h0+1e-6]=(*it);
			}
			else
				height2node[h0]=(*it);
		}
    }
	//map <Node*,int> node2rank;//del-loc
	//map <int,Node*> rank2node;//del-loc
	int rank=0;
	
	map <Node*,scalar_type> new_height;//del-loc
	int n=S->getNumberOfLeaves();
	scalar_type rank_height=0;
	
	for (map <scalar_type,Node *>::iterator hit=height2node.begin();hit!=height2node.end();hit++)
    {
		rank_height+=1.0/(scalar_type)(n-rank);
		rank+=1;
		//rank2node[rank]=(*hit).second;
		//node2rank[(*hit).second]=rank;
		new_height[(*hit).second]=rank_height;
    }
	
	for (vector <Node *>::iterator it=nodes.begin();it!=nodes.end();it++)
		if ((*it)->hasFather())
		{
			if ((*it)->isLeaf())
			{
				scalar_type fathers_height= new_height[ (*it)->getFather() ];
				(*it)->setDistanceToFather(fathers_height/rank_height);
			}
			else
			{
				scalar_type fathers_height= new_height[ (*it)->getFather() ];
				(*it)->setDistanceToFather((fathers_height- new_height[ (*it)])/rank_height);
			}
		}
	node2height.clear();
	height2node.clear();
	
	//node2rank.clear();
	//rank2node.clear();
	
	new_height.clear();
	
	nodes.clear();
	
	return ;
}


