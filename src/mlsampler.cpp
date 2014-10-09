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


  model->set_model_parameter("min_D",3);
  model->set_model_parameter("grid_delta_t",0.05);

  model->construct(Sstring);
  model->set_model_parameter("event_node",0);
  model->set_model_parameter("leaf_events",1);
  model->set_model_parameter("N",1);

  model->construct(Sstring);

  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("sigma_hat", 1);
  model->set_model_parameter("leaf_events",1);

  model->calculate_EGb();
  scalar_type old_p=model->p(ale);

  boost::timer * t = new boost::timer();
  cout <<"LL: " << log(model->p(ale)) << endl;
  cout <<"time: " << t->elapsed() << endl;
  cout << ".."<<endl; 
  vector<Tree*> sample_trees;
  string outname=ale_file+".rate_sample.samples";
  ofstream fout( outname.c_str() );
  string outname2=ale_file+".rate_sample.Ttokens";
  ofstream fout2( outname2.c_str() );

  int subsamples=atoi(argv[3]);
  boost::progress_display pd( subsamples );

  for (int i=0;i<90;i++) 
    {
      vector<scalar_type> ds;
      //rate proposal
      for (int i=0;i<3;i++)
	{
	  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
	  scalar_type d;
	  if (r<1./3.) d=RandomTools::randExponential(0.001)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
	  else if (r<2./3.) d=RandomTools::randExponential(0.01)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
	  else d=RandomTools::randExponential(0.1)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
	  ds.push_back(d);
	}
      scalar_type new_delta=delta+ds[0];
      scalar_type new_tau=tau+ds[1];
      scalar_type new_lambda=lambda+ds[2];
      
      //boundaries
      if (new_delta<1e-6) new_delta=1e-6;
      if (new_delta>10-1e-6) new_delta=10-1e-6;
      if (new_tau<1e-6) new_tau=1e-6;
      if (new_tau>10-1e-6) new_tau=10-1e-6;
      if (new_lambda<1e-6) new_lambda=1e-6;
      if (new_lambda>10-1e-6) new_lambda=10-1e-6;

      //likelihood
      model->set_model_parameter("delta",new_delta);
      model->set_model_parameter("tau",new_tau);
      model->set_model_parameter("lambda",new_lambda);
      model->calculate_EGb();
      scalar_type new_p=model->p(ale);
      if (new_p>=old_p or new_p/old_p>RandomTools::giveRandomNumberBetweenZeroAndEntry(1))
	{
	  old_p=new_p;
	  delta=new_delta; tau=new_tau; lambda=new_lambda;
	}
      cout << delta << " " << tau << " " << lambda << " " << log(old_p) << endl;

    }
  for (int i=0;i<subsamples;i++) 
    {		

      for (int i=0;i<10;i++) 
	{
	  vector<scalar_type> ds;
	  //rate proposal
	  for (int i=0;i<3;i++)
	    {
	      scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
	      scalar_type d;
	      if (r<1./3.) d=RandomTools::randExponential(0.001)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
	      else if (r<2./3.) d=RandomTools::randExponential(0.01)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
	      else d=RandomTools::randExponential(0.1)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
	      ds.push_back(d);
	    }
	  scalar_type new_delta=delta+ds[0];
	  scalar_type new_tau=tau+ds[1];
	  scalar_type new_lambda=lambda+ds[2];
      
	  //boundaries
	  if (new_delta<1e-6) new_delta=1e-6;
	  if (new_delta>10-1e-6) new_delta=10-1e-6;
	  if (new_tau<1e-6) new_tau=1e-6;
	  if (new_tau>10-1e-6) new_tau=10-1e-6;
	  if (new_lambda<1e-6) new_lambda=1e-6;
	  if (new_lambda>10-1e-6) new_lambda=10-1e-6;

	  //likelihood
	  model->set_model_parameter("delta",new_delta);
	  model->set_model_parameter("tau",new_tau);
	  model->set_model_parameter("lambda",new_lambda);
	  model->calculate_EGb();
	  scalar_type new_p=model->p(ale);
	  if (new_p>=old_p or new_p/old_p>RandomTools::giveRandomNumberBetweenZeroAndEntry(1))
	    {
	      old_p=new_p;
	      delta=new_delta; tau=new_tau; lambda=new_lambda;
	    }
	  cout << delta << " " << tau << " " << lambda << " " << log(old_p) << endl;
	}
      ++pd;
      string sample_tree=model->sample(false);
      fout << sample_tree << endl;
      for (vector<string>::iterator it=model->Ttokens.begin();it!=model->Ttokens.end();it++) fout << i << " " <<(*it)<<endl;

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

  string voutname=ale_file+".rate_sample.vstrings"; 
  ofstream vout( voutname.c_str() );

  for (std::map<long int, std::vector<int> >::iterator it=model->gid_branches.begin();it!=model->gid_branches.end();it++)
    {
      long int g_id=(*it).first;
      vout << g_id << " " << model->vertical_string(g_id) << endl;

    };

  return 1;
  

}

