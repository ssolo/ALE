#include "ALE.h"

using namespace std;
using namespace bpp;
#include <boost/algorithm/string.hpp>
#include <numeric>


boost::dynamic_bitset<> bpname(Node * node)
{
  Node * root=node;
  while (root->hasFather()) root=root->getFather();			      

  vector<string> all_leaves = TreeTemplateTools::getLeavesNames(*root);		
  vector<string> node_leaves = TreeTemplateTools::getLeavesNames(*node);
  set<string> node_set(node_leaves.begin(),node_leaves.end()); 
  sort(all_leaves.begin(), all_leaves.end());

  boost::dynamic_bitset<> bpname;
  for (int i=0;i<all_leaves.size();i++)
    {
      if (node_set.count(all_leaves[i]))
	bpname.push_back(true);
      else
	bpname.push_back(false);       	
    }
  return bpname;
}

string sorted_name(vector<string> left,vector<string> right)
{
  sort(left.begin(), left.end());
  sort(right.begin(), right.end());
  
  string sorted_name;
  if (left[0] > right[0])    
    sorted_name=left[0]+"\t"+right[0];
  else
    sorted_name=right[0]+"\t"+left[0];
  
  return sorted_name;
}

string azname(Node * node)
{
  if (node->isLeaf()) return node->getName();
    
  vector<string> left_aznames = TreeTemplateTools::getLeavesNames(*(node->getSons()[0]));
  vector<string> right_aznames = TreeTemplateTools::getLeavesNames(*(node->getSons()[1]));

  sort(left_aznames.begin(), left_aznames.end());
  sort(right_aznames.begin(), right_aznames.end());
  
  string azname=sorted_name(left_aznames,right_aznames);
  
  return azname;
}


int main(int argc, char ** argv)
{
  string line;
  ifstream tree_stream (argv[1]);

  map < boost::dynamic_bitset<> , Node * > RT_nodes;
    
  ifstream rooted_tree_stream (argv[2]);
  string tree;
  getline (rooted_tree_stream,tree);
  tree_type * RT=TreeTemplateTools::parenthesisToTree(tree,true);
  vector <Node *> nodes=RT->getNodes();
  for (auto it=nodes.begin();it!= nodes.end();it++)
    if ((*it)->hasFather())
      {	
	boost::dynamic_bitset<> bp=bpname((*it));
	//cout << bp << " " << azname((*it)) << endl;	
	RT_nodes[bp]=(*it);
      }
  map <string,vector <scalar_type> > bls;
  map <string, Node *> namenodes;
  string last_tree;
  tree_type * T;

  while (! tree_stream.eof())
    {
      getline (tree_stream,tree);

      if (tree.find("(")!=tree.npos) 
      {
	  vector <string> tokens;
	  boost::split(tokens,tree,boost::is_any_of(" \t"),boost::token_compress_on);
	  //cout <<"ur: "<< tokens[tokens.size()-1] << endl;
	  T=TreeTemplateTools::parenthesisToTree(tokens[tokens.size()-1],false);
	  vector <Node *> nodes=T->getNodes();
	  
	  map <scalar_type, string > names;           
	  for (auto it=nodes.begin();it!=nodes.end();it++)
	    if ((*it)->hasFather())
	      {
		Node * node = (*it);
		string name=azname(node);
		namenodes[name]=node;
		scalar_type bl=node->getDistanceToFather();
		boost::dynamic_bitset<> bp=bpname((*it));	  

		if (RT_nodes.count(bp)==1)
		  RT_nodes[bp]->setDistanceToFather(bl);
		if (RT_nodes.count(~bp)==1)
		  RT_nodes[~bp]->setDistanceToFather(bl);
		
		  
		if (bls.count(name)==0)
		  {
		    vector <scalar_type> tmp;
		    bls[name]=tmp;
		  }
		bls[name].push_back(bl);	       
	      }
	  //cout <<"RT: "<< TreeTemplateTools::treeToParenthesis(*RT,true);       	
      }

    }
  map <Node*, scalar_type> node_means;
  map <Node*, scalar_type> node_vars;
  
  for (map <string,vector <scalar_type> >::iterator nit=bls.begin();nit!=bls.end();nit++)
    {
      string name= (*nit).first ;
      scalar_type sum=0;
      scalar_type sq_sum=0;
      scalar_type n=0;
      
      for (int i=0;i<(*nit).second.size();i++)
	{
	  sum+=(*nit).second[i];
	  sq_sum+= (*nit).second[i] * (*nit).second[i];
	  n+=1;
	}
      scalar_type mean = sum/n;
      scalar_type var  = (sq_sum-sum*sum/n)/(n-1);
      node_means[namenodes[name]]=mean;
      node_vars[namenodes[name]]=var;
      
      //cout << "\t"<< (*nit).second[i];
      //cout << endl;
    }

  nodes=T->getNodes();
  for (auto it=nodes.begin();it!=nodes.end();it++)
    if ((*it)->hasFather())
      {
	Node * node = (*it);
	//node->setDistanceToFather(node_means[node]);
	boost::dynamic_bitset<> bp=bpname(node);	  
	
	if (RT_nodes.count(bp)==1)
	  RT_nodes[bp]->setDistanceToFather(node_means[node]);
	if (RT_nodes.count(~bp)==1)
	  RT_nodes[~bp]->setDistanceToFather(node_means[node]);
	
      }
  cout << TreeTemplateTools::treeToParenthesis(*RT,true);
  for (auto it=nodes.begin();it!=nodes.end();it++)
    if ((*it)->hasFather())
      {
	Node * node = (*it);
	//node->setDistanceToFather(node_vars[node]);
	boost::dynamic_bitset<> bp=bpname(node);	  
	
	if (RT_nodes.count(bp)==1)
	  RT_nodes[bp]->setDistanceToFather(node_vars[node]*1e6);
	if (RT_nodes.count(~bp)==1)
	  RT_nodes[~bp]->setDistanceToFather(node_vars[node]*1e6);
	
      }
  cout << TreeTemplateTools::treeToParenthesis(*RT,true);

      
  return 1;
  
}
