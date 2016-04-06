#include "ALE_util.h"
#include "mpi_tree.h"

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>

using namespace std;
using namespace bpp;
using namespace boost::mpi;

string random_tree_newick(string Sstring)
{
  tree_type * T=TreeTemplateTools::parenthesisToTree(Sstring,false);
  vector<string> random_tree_population; 
  map<string,scalar_type> random_tree_ages;
  vector <string> leaf_names=T->getLeavesNames();
    for (vector<string>::iterator it=leaf_names.begin();it!=leaf_names.end();it++) 
      {
	random_tree_population.push_back((*it));
	random_tree_ages[(*it)]=0;
      }
  scalar_type t=0;
  while(random_tree_population.size()>1)
    {
      int Nr=random_tree_population.size();
      scalar_type t_next=RandomTools::randExponential(1./(2*Nr));
      t+=t_next;
      int i=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(Nr);    
      int j=i;
      while (i==j) j=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(Nr);
      stringstream tmp;
      tmp<<"("<<random_tree_population[i]<<":"<<t-random_tree_ages[random_tree_population[i]]<<","<<random_tree_population[j]<<":"<<t-random_tree_ages[random_tree_population[j]]<<")";       
      random_tree_ages[tmp.str()]=t;
      if (j>i)
	{
	  random_tree_population.erase(random_tree_population.begin()+j);
	  random_tree_population.erase(random_tree_population.begin()+i);
	}
      else
	{
	  random_tree_population.erase(random_tree_population.begin()+i);
	  random_tree_population.erase(random_tree_population.begin()+j);
	}
      random_tree_population.push_back(tmp.str());	 
    }

  //cout << random_tree_population[0]<<";" << endl;

  return random_tree_population[0]+";"; 
}


int main(int argc, char ** argv)
{
  map <string,scalar_type > ll_cache;

  environment env(argc, argv);
  communicator world;
  int done=1;
  int it_num=2;
  bool bw=1;
  int N_SPR=50;
  ifstream file_stream_S (argv[1]);

  string Sstring;
  getline (file_stream_S,Sstring);
  string Rstring=Sstring;
  if (atoi(argv[3])==1) Rstring=random_tree_newick(Sstring);
  map<string,scalar_type> parameters;
  mpi_tree * infer_tree = new mpi_tree(Sstring,world,parameters,true);
  infer_tree->load_distributed_ales(argv[2]);

  scalar_type Sll = infer_tree->calculate_pun(it_num);
	    

  string old_S=Rstring;
  string max_S=Rstring;
  infer_tree->model->construct_undated(max_S);
  scalar_type max_ll = infer_tree->calculate_pun(it_num,bw);
  scalar_type new_ll;
  infer_tree->gather_T_to_from();

  //infer_tree->model->construct_undated(max_S);
  //infer_tree->calculate_pun(it_num,1);
  //infer_tree->gather_T_to_from();
  
  bool changed=true;
  bool last=false;
  while (changed)
    {
      changed=false;
      infer_tree->model->construct_undated(max_S);

      if (world.rank()==0) cout <<"@ll "<< max_ll << " " << Sll << " " << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(max_S),*TreeTemplateTools::parenthesisToTree(Sstring)) << endl;


      string old_S=max_S;
      
      vector <int> sorted_e;
      vector <int> sorted_f;
      vector <int> sorted_Ts;

      if (world.rank()==0)
	for (map <scalar_type, vector< int > >::iterator it=infer_tree->sort_e.begin();it!=infer_tree->sort_e.end();it++)
	  {
	    scalar_type Ts=(*it).first;
	    for (int i=0;i<(*it).second.size();i++)
	      {
		int e=infer_tree->sort_e[Ts][i];
		int f=infer_tree->sort_f[Ts][i];
		sorted_e.push_back(e);
		sorted_f.push_back(f);
		sorted_Ts.push_back(Ts);				
	      }
	  }
      broadcast(world,sorted_e,0);
      broadcast(world,sorted_f,0);
      for (int i=0;i<N_SPR;i++)
	{
	  infer_tree->model->construct_undated(old_S);

	  string new_S=infer_tree->model->feSPR(sorted_e[i],sorted_f[i]);
	  if (ll_cache.count(new_S)==0)
	    {
	      infer_tree->model->construct_undated(new_S);	  
	      new_ll= infer_tree->calculate_pun(it_num,bw);
	      ll_cache[new_S]=new_ll;
	    }
	  else
	    {
	      new_ll=ll_cache[new_S];
	    }
	  if (world.rank()==0)
	    {
	      if (sorted_e[i]<infer_tree->model->last_leaf)
		cout << " " << infer_tree->model->node_name[infer_tree->model->id_nodes[sorted_e[i]]];
	      else
		cout << " " << sorted_e[i];
	      if (sorted_f[i]<infer_tree->model->last_leaf)
		cout << "->" << infer_tree->model->node_name[infer_tree->model->id_nodes[sorted_f[i]]];	      
	      else
		cout << "->" << sorted_f[i];
	      cout << " with: " << -sorted_Ts[i]<< " Ts " ; //" " << new_S << endl;
	    }
 	  if (world.rank()==0) cout << new_ll << " " << max_ll << " " << Sll << " " 
				    << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(Sstring)) <<" "
				    << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(Rstring)) << " "
				    << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(old_S)) << endl;	  
	  if ( (new_ll>max_ll) )
	    {
	      infer_tree->gather_T_to_from();
	      max_S=new_S;
	      changed=true;
	      max_ll=new_ll;
	      break;
	    }
	}
      if (not changed)
	{
	  N_SPR=0;
	  infer_tree->model->construct_undated(old_S);
	  if (world.rank()==0) cout<< " .. new roots .. "  << infer_tree->model->string_parameter["S_with_ranks"] << endl;
	  int e=infer_tree->model->last_branch-1;
	  {
	    infer_tree->model->construct_undated(old_S);      
	    vector<string> NNIs=infer_tree->model->NNIs(e);
	    for (int nni=0;nni<NNIs.size();nni++)
	      {	      
		string new_S=NNIs[nni];
		if (ll_cache.count(new_S)==0)
		  {
		    infer_tree->model->construct_undated(new_S);	  
		    new_ll= infer_tree->calculate_pun(it_num,bw);
		    ll_cache[new_S]=new_ll;
		  }
		else
		  {
		    new_ll=ll_cache[new_S];
		  }
		if (world.rank()==0) cout << new_ll << " " << max_ll << " " << Sll << " " 
					  << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(Sstring)) <<" "
					  << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(Rstring)) << " "
					  << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(old_S)) << endl;	       
		if ( (new_ll>max_ll) )
		    {
		      max_S=new_S;
		      max_ll=new_ll;
		      changed=true;
		      break;
		    }
		}
	    if (changed) e=0;
	  }
	}
      if (not changed)
	{
	  infer_tree->model->construct_undated(old_S);
	  if (world.rank()==0) cout<< " .. all NNIs .. "  << infer_tree->model->string_parameter["S_with_ranks"] << endl;
	  for (int e=infer_tree->model->last_branch-2;e>infer_tree->model->last_leaf-1;e--)
	    {
	      infer_tree->model->construct_undated(old_S);      
	      vector<string> NNIs=infer_tree->model->NNIs(e);
	      for (int nni=0;nni<NNIs.size();nni++)
		{	      
		  string new_S=NNIs[nni];
		  if (ll_cache.count(new_S)==0)
		    {
		      infer_tree->model->construct_undated(new_S);	  
		      new_ll= infer_tree->calculate_pun(it_num,bw);
		      ll_cache[new_S]=new_ll;
		    }
		  else
		    {
		      new_ll=ll_cache[new_S];
		    }
		  if (world.rank()==0) cout << new_ll << " " << max_ll << " " << Sll << " " 
					    << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(Sstring)) <<" "
					    << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(Rstring)) << " "
					    << TreeTools::robinsonFouldsDistance(*TreeTemplateTools::parenthesisToTree(new_S),*TreeTemplateTools::parenthesisToTree(old_S)) << endl;	  		  
		  if ( (new_ll>max_ll) )
		    {
		      max_S=new_S;
		      max_ll=new_ll;
		      changed=true;
		      break;
		    }
		}	      
	      if (changed) e=0;
	    }
	}
      if (world.rank()==0) cout<< "@ "  << infer_tree->model->string_parameter["S_with_ranks"] << endl;
      if (not changed and not last)
	{
	  last=true;
	  changed=true;
	  N_SPR=200;
	}
      broadcast(world,done,0);	      
    }
  if (world.rank()==0) cout<< ". "  << infer_tree->model->string_parameter["S_with_ranks"] << endl;

}

