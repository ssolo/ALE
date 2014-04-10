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

  string ml_rec_file=argv[2];

  vector <string> tokens;
  boost::split(tokens,ml_rec_file,boost::is_any_of("."),boost::token_compress_on);

  string ale_file=tokens[0]+".ale";
  cout << ale_file << endl;

  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);

  scalar_type delta,tau,lambda;
  ifstream ml_file_stream (ml_rec_file.c_str());
  while(! ml_file_stream.eof())
    {
      string line;
      getline (ml_file_stream,line);
      if (line.find("ML")!=line.npos )
	{
	  vector <string> tokens;
	  boost::split(tokens,line,boost::is_any_of(" \t"),boost::token_compress_on);
	  delta=atof(tokens[1].c_str());
	  tau=atof(tokens[2].c_str());
	  lambda=atof(tokens[3].c_str());

	}
    }
  cout << delta << " " << tau << " " << lambda << endl;
  exODT_model* model=new exODT_model();


  model->set_model_parameter("min_D",4);
  model->set_model_parameter("grid_delta_t",0.005);
  model->set_model_parameter("DD",10);

  model->construct(Sstring);

  model->set_model_parameter("event_node",0);  
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("leaf_events",1);

  model->calculate_EGb();

  boost::timer * t = new boost::timer();
  cout <<"LL: " << model->p(ale) << endl;
  cout <<"time: " << t->elapsed() << endl;
  cout << ".."<<endl; 
  vector<Tree*> sample_trees;
  string outname=ale_file+".samples";
  ofstream fout( outname.c_str() );
  string outname2=ale_file+".Ttokens";
  ofstream fout2( outname2.c_str() );

  int subsamples=atoi(argv[3]);
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

