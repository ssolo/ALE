#include "exODT.h"
#include "exODT_sim.h"

#include "ALE_util.h"
#include <boost/progress.hpp>

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  //we need a species tree

  string ml_rec_file=argv[1];
  string Sstring;
  vector <string> tokens;
  boost::split(tokens,ml_rec_file,boost::is_any_of("."),boost::token_compress_on);
  string ale_file=argv[2];

  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);

  scalar_type delta,tau,lambda;
  ifstream ml_file_stream (ml_rec_file.c_str());
  while(! ml_file_stream.eof())
    {
      string line;
      getline (ml_file_stream,line);
      if (line.find("S:")!=line.npos )
	{
	  vector <string> tokens;
	  boost::split(tokens,line,boost::is_any_of(" \t"),boost::token_compress_on);
	  Sstring=tokens[1];
	}
      if (line.find("ML")!=line.npos )
	{
	  vector <string> tokens;
	  boost::split(tokens,line,boost::is_any_of(" \t"),boost::token_compress_on);
	  delta=atof(tokens[1].c_str());
	  tau=atof(tokens[2].c_str());
	  lambda=atof(tokens[3].c_str());

	}
    }
  cout <<"# S: "<< Sstring << endl;
  cout <<"# rates: "<< delta << " " << tau << " " << lambda << endl;
  exODT_model* model=new exODT_model();

  model->set_model_parameter("BOOTSTRAP_LABELS","yes");

  model->set_model_parameter("min_D",3);
  model->set_model_parameter("grid_delta_t",0.05);
  model->construct(Sstring);

  model->set_model_parameter("event_node",0);
  model->set_model_parameter("leaf_events",1);
  model->set_model_parameter("N",1);


  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("sigma_hat", 1);

  cout << "." << endl; 
  model->calculate_EGb();
  cout << ".." << endl; 

  boost::timer * t = new boost::timer();
  cout <<"LL: " << log(model->p(ale)) << endl;
  cout <<"time: " << t->elapsed() << endl;
  cout << ".."<<endl; 
  vector<Tree*> sample_trees;
<<<<<<< HEAD
  string outname=ml_rec_file+".rate_resample.samples";
  ofstream fout( outname.c_str() );
  string outname2=ml_rec_file+".rate_resample.Ttokens";
=======
  string outname=ale_file+".rate_resample.samples";
  ofstream fout( outname.c_str() );
  string outname2=ale_file+".rate_resample.Ttokens";
>>>>>>> 474d82cc0794533095baad28e4c97d06f7545120
  ofstream fout2( outname2.c_str() );

  int samples=atoi(argv[3]);
  //boost::progress_display pd( samples );
  for (int i=0;i<samples;i++)
    {
      //++pd;
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

  scalar_type cumsum=0;
  map<scalar_type,scalar_type> cum_copies;
  for (int branch=model->last_branch-1;branch>=0;branch--)	
    {
      scalar_type extant;
      if (model->id_ranks[branch]==0)
	extant=model->last_branch;
      else
	extant=model->id_ranks[branch]-model->last_branch;
      cumsum+=model->branch_counts["copies"][branch];
      //cum_copies[model->branch_ts[branch]]=cumsum;
      cum_copies[model->id_ranks[branch]]=cumsum;
    }
  for (map<scalar_type,scalar_type>::iterator it=cum_copies.begin();it!=cum_copies.end();it++)
    cout << (*it).first << " " << (*it).second << endl;

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
  return 1;
  pair<string, scalar_type> res = model->p_MLRec(ale);    
  cout << endl;
  cout << "ML: "<< endl; 
  cout << res.first << endl;
  cout << endl;

<<<<<<< HEAD
  string voutname=ml_rec_file+".rate_resample.vstrings"; 
=======
  string voutname=ale_file+".rate_resample.vstrings"; 
>>>>>>> 474d82cc0794533095baad28e4c97d06f7545120
  ofstream vout( voutname.c_str() );

  for (std::map<long int, std::vector<int> >::iterator it=model->gid_branches.begin();it!=model->gid_branches.end();it++)
    {
      long int g_id=(*it).first;
      vout << g_id << " " << model->vertical_string(g_id) << endl;

    };

  return 1;
  

}

