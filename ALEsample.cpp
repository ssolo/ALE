
#include "exODT.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;

int main(int argc, char ** argv)
{
  //RandomTools::setSeed(20110426);
  
  //we need a species tree
  string Sstring;
  ifstream file_stream (argv[1]);
  getline (file_stream,Sstring);

  //we need an .ale file containing observed conditional clade probabilities
  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);

  //we initialise a coarse grained reconciliation model for calculating the sum
  exODT_model* model=new exODT_model();
  //
  model->set_model_parameter("D",1);
  model->set_model_parameter("DD",10);
  model->construct(Sstring);
  //cf. section S1.4 of the Supporting Material of www.pnas.org/cgi/doi/10.1073/pnas.1202997109
  //for a description of the event node approximation 
  model->set_model_parameter("event_node",1);

  //and a finer grained reconciliation model for sampling
  sample_model->set_model_parameter("D",10);
  sample_model->set_model_parameter("DD",10);
  sample_model->construct(Sstring);
  //stocahstic sampling is not compatible with the event node approximation
  sample_model->set_model_parameter("event_node",0);
  
  sample_model->set_model_parameter("leaf_events",1);

  //a set of inital rates
  scalar_type delta=0.01,tau=0.01,lambda=0.02;
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau",tau);
  model->set_model_parameter("lambda",lambda);
  //calculate_EGb() must always be called after chaging rates to calculate E-s and G-s 
  //cf. http://arxiv.org/abs/1211.4606
  model->calculate_EGb();

  //p(ale) calculates the probability of the obsorved set of gene trees, i.e. Pi(Gamma) 
  //cf. ALEPAPER
  scalar_type old_p=model->p(ale);

  int steps=0;
  int N_steps=1000;
  int burnin=100;
  bool allprint=false;
  int print_mode=1;
  int samples=10;
  vector<Tree*> sample_trees;
  vector<string> sample_strings;


  while (steps<N_steps+burnin)
    {
      vector<scalar_type> ds;
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
      if (new_delta<1e-6) new_delta=1e-6;
      if (new_delta>10-1e-6) new_delta=10-1e-6;
      if (new_tau<1e-6) new_tau=1e-6;
      if (new_tau>10-1e-6) new_tau=10-1e-6;
      if (new_lambda<1e-6) new_lambda=1e-6;
      if (new_lambda>10-1e-6) new_lambda=10-1e-6;
      model->set_model_parameter("delta",new_delta);
      model->set_model_parameter("tau",new_tau);
      model->set_model_parameter("lambda",new_lambda);
      model->calculate_EGb();
      scalar_type new_p=model->p(ale);
      //cout << "#try " <<  new_delta << " " << new_tau << " " << new_lambda << " " << new_p/old_p << endl; 
      if (new_p>=old_p or new_p/old_p>RandomTools::giveRandomNumberBetweenZeroAndEntry(1))
	{
	  old_p=new_p;
	  delta=new_delta; tau=new_tau; lambda=new_lambda;
	  steps++;
	  if ( ( steps>burnin and steps%print_mode==0 ) or allprint) 
	    {
	      cout << "sample_rates "<< steps << "\t" << delta << "\t" << tau << "\t" << lambda << "\t" << log(old_p) << endl; 
	      for (int i=0;i<samples;i++) 
		{
		  string sample_tree=model->sample(false);
		  sample_trees.push_back(TreeTemplateTools::parenthesisToTree(sample_tree));
		  sample_strings.push_back(sample_tree);
		  cout << "sample_rec " << steps << "\t" << i << "\t" << sample_tree << endl;
		  cout << "sample_events " << steps << "\t" << i << "\t" << model->MLRec_events["D"] << "\t" << model->MLRec_events["T"] << "\t" << model->MLRec_events["L"]<< "\t" << model->MLRec_events["S"] <<endl;
		}	      
	    }
	}      
    }
  cout << endl;
  Tree* con_tree= TreeTools::thresholdConsensus(sample_trees,0.5);
  cout << TreeTemplateTools::treeToParenthesis(*con_tree) << endl;
  TreeTools::computeBootstrapValues(*con_tree,sample_trees);
  string con_tree_sup=TreeTemplateTools::treeToParenthesis(*con_tree);
  cout << "#cons. tree: " << con_tree_sup;

}
