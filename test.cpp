#include "exODT.h"
#include "exODT_sim.h"

#include "ALE_util.h"
#include <omp.h>
#include <boost/progress.hpp>

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  //we need a species tree

  string sname=argv[1];

  string Sstring;
  ifstream file_stream (sname.c_str());
  getline (file_stream,Sstring);

  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);

  exODT_model* model=new exODT_model();

  //XX
  //XX

  //exODT_sim* simulation=new exODT_sim(100,1010);
  
  scalar_type delta=atof(argv[3]);
  scalar_type tau=atof(argv[4]);
  scalar_type lambda=atof(argv[5]);

  //simulation->sample_species(10);

  
  //for (vector<string>::iterator it=simulation->gene_trees.begin();it!=simulation->gene_trees.end();it++)
  //cout << (*it) << endl;

  //cout << simulation->S_string << endl;

  model->set_model_parameter("min_D",3);
  model->set_model_parameter("grid_delta_t",0.005);
  model->set_model_parameter("DD",10);

  model->construct(Sstring);

  model->set_model_parameter("event_node",0);  
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("leaf_events",1);

  model->calculate_EGb();
  cout << model->p(ale) << endl;
  cout << ".."<<endl; 

  /*
  pair<string, scalar_type> res = model->p_MLRec(ale);    
  cout << res.first <<endl;
  cout << endl;
  cout << res.second <<endl;
  cout << endl;
  cout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
  cout <<"Total \t"<< model->MLRec_events["D"] << "\t" << model->MLRec_events["T"] << "\t" << model->MLRec_events["L"]<< "\t" << model->MLRec_events["S"] <<endl;     
  */
  vector<Tree*> sample_trees;
  string outname=ale_file+".samples";
  ofstream fout( outname.c_str() );
  string outname2=ale_file+".Ttokens";
  ofstream fout2( outname2.c_str() );

  int subsamples=atoi(argv[6]);
  boost::progress_display pd( subsamples );

  for (int i=0;i<subsamples;i++) 
    {		  
      ++pd;
      string sample_tree=model->sample(false);
      fout << sample_tree << endl;
      for (vector<string>::iterator it=model->Ttokens.begin();it!=model->Ttokens.end();it++) fout2<< i << " " <<(*it)<<endl;

      tree_type * G=TreeTemplateTools::parenthesisToTree(sample_tree,false);
      vector<Node*> leaves = G->getLeaves();
      for (vector<Node*>::iterator it=leaves.begin();it!=leaves.end();it++ )
	{
	  string name=(*it)->getName();
	  vector<string> tokens;
	  boost::split(tokens,name,boost::is_any_of(".@"),boost::token_compress_on);
	  (*it)->setName(tokens[0]);
	  tokens.clear();
	}
      leaves.clear();
      sample_trees.push_back(G);	      
    }
  
  cout << model->counts_string();
  cout << "Os" <<endl;
  model->show_counts("Os");
  cout << "Ds" <<endl;
  model->show_counts("Ds");
  cout << "Ts" <<endl;
  model->show_counts("Ts");
  cout << "Ts from" <<endl;
  model->show_counts("Tfroms");
  cout << "Ls" <<endl;
  model->show_counts("Ls");
  cout << "copies" <<endl;
  model->show_counts("copies");

  Tree* con_tree= TreeTools::thresholdConsensus(sample_trees,0.5);
  TreeTools::computeBootstrapValues(*con_tree,sample_trees);
  cout << endl;
  cout << "thcon: "<<endl;
  string con_str = TreeTemplateTools::treeToParenthesis(*con_tree);    
  cout << con_str;
  cout << endl;


  /*  
  approx_posterior * sale=observe_ALE_from_strings(sample_strings);
  pair<string,scalar_type> mpp_res=sale->mpp_tree();
  Tree* mpp_T = TreeTemplateTools::parenthesisToTree(mpp_res.first,false);
  TreeTools::computeBootstrapValues(*mpp_T,sample_trees);
  cout << endl;
  cout << "mpp: "<<endl;
  cout << mpp_res.first << endl;
  //cout << TreeTools::treeToParenthesis(*mpp_T);
  cout << endl;
  approx_posterior * cale=observe_ALE_from_string(con_str);
  */

  pair<string, scalar_type> res = model->p_MLRec(ale);    
  cout << endl;
  cout << "ML: "<< endl; 
  cout << res.first << endl;
  cout << endl;

  
  return 1;
  

}

