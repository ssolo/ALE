
#include "exODT.h"
#include "ALE_util.h"
#include <boost/progress.hpp>

using namespace std;
using namespace bpp;

int main(int argc, char ** argv)
{
  cout << "ALEsample using ALE v"<< ALE_VERSION <<endl;

  //RandomTools::setSeed(20110426);
  if (argc<3) 
    {
      cout << "usage:\n ./ALEsample species_tree.newick gene_tree_sample.ale [samples=1000] [burnin=100]" << endl;
      return 1;
    }
  
  //we need a dared species tree in newick format
  string Sstring;
  ifstream file_stream_S (argv[1]);
  getline (file_stream_S,Sstring);
  cout << "Read species tree from: " << argv[1] <<".."<<endl;
  //we need an .ale file containing observed conditional clade probabilities
  //cf. ALEobserve
  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);
  cout << "Read summary of tree sample for "<<ale->observations<<" trees from: " << ale_file <<".."<<endl;

  //we initialise the model 
  exODT_model* model=new exODT_model();
  //
  model->set_model_parameter("BOOT_STRAP_LABLES","yes");

  model->set_model_parameter("min_D",3);
  model->set_model_parameter("grid_delta_t",0.05);
  model->construct(Sstring);

  model->set_model_parameter("event_node",0);
  model->set_model_parameter("leaf_events",1);
  model->set_model_parameter("N",1);

  //a set of inital rates
  scalar_type delta=0.01,tau=0.01,lambda=0.1;
  string append="";
  if (argc>4)
    append=argv[4];
  if (argc>7)
    delta=atof(argv[5]),tau=atof(argv[6]),lambda=atof(argv[7]);  
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("sigma_hat", 1);

  model->calculate_EGb();

  //p(ale) calculates the probability of the obsorved set of gene trees, i.e. Pi(Gamma) 
  //cf. ALEPAPER
  scalar_type old_p=model->p(ale);

  int steps=0,accapted=0,sampled=0;
  int N_samples=1000;
  if (argc>3) N_samples=atoi(argv[3]);
  int burnin=100;
  if (argc>4) burnin=atoi(argv[4]);
  bool allprint=false;
  int print_mod=10;
  int subsamples=10;
  vector<Tree*> sample_trees;
  vector<string> sample_strings;

  cout << "Starting burnin.."  << endl;

  string rate_name=ale_file+append+".ratelist";
  string event_name=ale_file+append+".eventlist";
  string sample_name=ale_file+append+".treelist";

  ofstream rate_out( rate_name.c_str() );
  ofstream event_out( event_name.c_str() );
  ofstream sample_out( sample_name.c_str() );
  
  rate_out <<  "#ALEsample using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl;
  rate_out <<  "#" << "sample" << "\t" << "step" << "\t" << "duplication_rate" << "\t" << "transfer_rate" << "\t" << "loss_rate" << "\t" << "log_likelihood"<< endl;
  
  event_out <<  "#ALEsample using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl;
  event_out <<  "#" << "subsample"<<"\t"<< "sample" << "\t" << "step" << "\t" << "duplications" << "\t" << "transfers" << "\t" << "losses" << "\t" << "speciations"<< endl;

  sample_out <<  "#ALEsample using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl;

  boost::progress_display pd( burnin );
  
  while (sampled<N_samples)
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
	  if (accapted<burnin) ++pd;
	  if (accapted==burnin) 
	    { 
	      cout << "\nFinished burnin. \n Sampling:" <<endl;   
	      pd.restart(N_samples);
	    }
	  accapted++;
	}      
      steps++;

      if ( ( accapted>burnin and steps%print_mod==0 ) or allprint) 
	{
	  sampled++;
	  ++pd;
	  model->set_model_parameter("delta",delta);
	  model->set_model_parameter("tau",tau);
	  model->set_model_parameter("lambda",lambda);
	  model->calculate_EGb();
	  old_p= model->p(ale);
	      
	  rate_out << sampled << "\t" << steps << "\t" << delta << "\t" << tau << "\t" << lambda << "\t" << log(old_p) << endl;
	  
	  for (int i=0;i<subsamples;i++) 
	    {		  
	      string sample_tree=model->sample(false);
	      sample_out << sample_tree << endl;
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
	      event_out << i << "\t" << sampled << "\t" << steps << "\t" << model->MLRec_events["D"] << "\t" << model->MLRec_events["T"] << "\t" << model->MLRec_events["L"]<< "\t" << model->MLRec_events["S"] <<endl;
	    }	      
	}      
    }
  cout << "Calculating consensus tree."<<endl;
  Tree* con_tree= TreeTools::thresholdConsensus(sample_trees,0.5);

  string con_name=ale_file+append+".cons_tree";

  ofstream con_out( con_name.c_str() );

  con_out <<  "#ALEsample using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl;
  TreeTools::computeBootstrapValues(*con_tree,sample_trees);
  string con_tree_sup=TreeTemplateTools::treeToParenthesis(*con_tree);
  cout <<endl<< "#cons. tree: " << con_tree_sup;

  cout << "Results in:" << endl;
  cout << rate_name << endl;
  cout << event_name << endl;
  cout << sample_name << endl;
  cout << con_name << endl;

}
