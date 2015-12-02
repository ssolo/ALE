
#include "exODT.h"
#include "ALE_util.h"

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>

using namespace std;
using namespace bpp;
 
class p_fun:
  public virtual Function,
  public AbstractParametrizable
{
private:
  double fval_;
  exODT_model* model_pointer;
  approx_posterior* ale_pointer;
public:
  p_fun(exODT_model* model,approx_posterior* ale, double delta_start=0.01,double tau_start=0.01,double lambda_start=0.1//,double sigma_hat_start=1.
) : AbstractParametrizable(""), fval_(0), model_pointer(model), ale_pointer(ale) 
  {
    //We declare parameters here:
 //   IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
      IntervalConstraint* constraint = new IntervalConstraint ( 1e-6, 10-1e-6, true, true );
      addParameter_( new Parameter("delta", delta_start, constraint) ) ;
      addParameter_( new Parameter("tau", tau_start, constraint) ) ;
      addParameter_( new Parameter("lambda", lambda_start, constraint) ) ;
      //addParameter_( new Parameter("sigma_hat", sigma_hat_start, constraint) ) ;

  }
  
  p_fun* clone() const { return new p_fun(*this); }
  
public:
  
    void setParameters(const ParameterList& pl)
    throw (ParameterNotFoundException, ConstraintException, Exception)
    {
        matchParametersValues(pl);
    }
    double getValue() const throw (Exception) { return fval_; }
    void fireParameterChanged(const ParameterList& pl)
    {
        double delta = getParameterValue("delta");
        double tau = getParameterValue("tau");
        double lambda = getParameterValue("lambda");
        //double sigma_hat = getParameterValue("sigma_hat");
        
        model_pointer->set_model_parameter("delta",delta);
        model_pointer->set_model_parameter("tau",tau);
        model_pointer->set_model_parameter("lambda",lambda);
        //model_pointer->set_model_parameter("sigma_hat",sigma_hat);
        model_pointer->calculate_undatedEs();
        double y=-log(model_pointer->pun(ale_pointer));
        //cout <<endl<< "delta=" << delta << "\t tau=" << tau << "\t lambda=" << lambda //<< "\t lambda="<<sigma_hat << "\t ll="
	//    << -y <<endl;
        fval_ = y;
    }
};


int main(int argc, char ** argv)
{
  cout << "ALEml using ALE v"<< ALE_VERSION <<endl;

  if (argc<3) 
    {
      cout << "usage:\n ./ALEml_undated species_tree.newick gene_tree_sample.ale [samples] [gene_name_seperator]" << endl;
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

  //we initialise a coarse grained reconciliation model for calculating the sum
  exODT_model* model=new exODT_model();

  scalar_type samples=100;
  if (argc>3)
    samples=atof(argv[3]);

  if (argc>4)
    model->set_model_parameter("gene_name_separators", argv[4]);
  model->set_model_parameter("BOOT_STRAP_LABLES","yes");

  model->construct_undated(Sstring);
  

  //a set of inital rates
  scalar_type delta=0.01,tau=0.01,lambda=0.1;  
  if (argc>7)
    delta=atof(argv[5]),tau=atof(argv[6]),lambda=atof(argv[7]);  
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);

  //calculate_EGb() must always be called after changing rates to calculate E-s and G-s
  //cf. http://arxiv.org/abs/1211.4606
  model->calculate_undatedEs();
  cout << "Reconciliation model initialised, starting DTL rate optimisation" <<".."<<endl;

  //we use the Nelder–Mead method implemented in Bio++
  Function* f = new p_fun(model,ale,delta,tau,lambda//,1
			  );
  Optimizer* optimizer = new DownhillSimplexMethod(f);

  optimizer->setProfiler(0);
  optimizer->setMessageHandler(0);
  optimizer->setVerbose(2);
  
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  optimizer->init(f->getParameters()); //Here we optimize all parameters, and start with the default values.
        
 
  optimizer->optimize();



  delta=optimizer->getParameterValue("delta");
  tau=optimizer->getParameterValue("tau");
  lambda=optimizer->getParameterValue("lambda");

  scalar_type mlll=-optimizer->getFunctionValue();
  cout << endl << "ML rates: " << " delta=" << delta << "; tau=" << tau << "; lambda="<<lambda//<<"; sigma="<<sigma_hat
       <<"."<<endl;
  
    

  cout << "LL=" << mlll << endl;
  
  cout << "Sampling reconciled gene trees.."<<endl;
  vector <string> sample_strings;
  vector <Tree*> sample_trees;
  boost::progress_display pd( samples );

  for (int i=0;i<samples;i++) 
    {
      ++pd;
      string sample_tree=model->sample_undated();
      sample_strings.push_back(sample_tree);      
      if (ale->last_leafset_id>3)
	{
	  
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
    }
  /*cout << "Calculating ML reconciled gene tree.."<<endl; 
  pair<string, scalar_type> res = model->p_MLRec(ale);    
  //and output it..
  */
  string outname=ale_file+".uml_rec"; 
  ofstream fout( outname.c_str() );
  fout <<  "#ALEml_undated using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl<<endl;
  fout << "S:\t"<<model->string_parameter["S_with_ranks"] <<endl;
  fout << endl;
  fout << "Input ale from:\t"<<ale_file<<endl;
  fout << ">logl: " << mlll << endl;
  fout << "rate of\t Duplications\tTransfers\tLosses" <<endl;
  fout << "ML \t"<< delta << "\t" << tau << "\t" << lambda //<< "'t" << sigma_hat
       << endl;
  fout << endl;
  fout << samples << " reconciled G-s:\n"<<endl;
  for (int i=0;i<samples;i++)
    {
      fout<<sample_strings[i]<<endl;
    }
  
  //fout << "reconciled G:\t"<< res.first <<endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
  fout <<"Total \t"<< model->MLRec_events["D"]/samples << "\t" << model->MLRec_events["T"]/samples << "\t" << model->MLRec_events["L"]/samples<< "\t" << model->MLRec_events["S"]/samples <<endl;    
  fout << endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tcopies" <<endl; 
  fout << model->counts_string_undated(samples);
  
  cout << "Results in: " << outname << endl;
  if (ale->last_leafset_id>3)
    {
      cout << "Calculating consensus tree."<<endl;
      Tree* con_tree= TreeTools::thresholdConsensus(sample_trees,0.5);
      
      string con_name=ale_file+".ucons_tree";
      
      ofstream con_out( con_name.c_str() );
      con_out <<  "#ALEsample using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl;
      TreeTools::computeBootstrapValues(*con_tree,sample_trees);
      string con_tree_sup=TreeTemplateTools::treeToParenthesis(*con_tree);
      con_out << con_tree_sup << endl;
      cout << endl<< "Consensus tree in " << con_name<< endl;
    }
  
  string t_name=ale_file+".uTs";
  ofstream tout( t_name.c_str() );
  tout <<"#from\tto\tfreq.\n";
  
  for (int e=0;e<model->last_branch;e++)
    for (int f=0;f<model->last_branch;f++)
      if  (model->T_to_from[e][f]>0)
	{
	  if (e<model->last_leaf)
	    tout << "\t" << model->node_name[model->id_nodes[e]];
	  else
	    tout << "\t" << e;
	  if (f<model->last_leaf)
	    tout << "\t" << model->node_name[model->id_nodes[f]];	      
	  else
	    tout << "\t" << f;
	  tout << "\t" << model->T_to_from[e][f]/samples <<  endl;
	}
  cout << "Transfers in: " << t_name << endl;
  return 0;
}
