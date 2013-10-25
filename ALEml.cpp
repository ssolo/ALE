
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
  p_fun(exODT_model* model,approx_posterior* ale, double delta_start=0.01,double tau_start=0.01,double lambda_start=0.1) : AbstractParametrizable(""), fval_(0), model_pointer(model), ale_pointer(ale) 
  {
    //We declare parameters here:
    IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
    addParameter_(Parameter("delta", delta_start,constraint));
    addParameter_(Parameter("tau", tau_start,constraint));
    addParameter_(Parameter("lambda", lambda_start,constraint));    
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

    model_pointer->set_model_parameter("delta",delta);
    model_pointer->set_model_parameter("tau",tau);
    model_pointer->set_model_parameter("lambda",lambda);
    model_pointer->calculate_EGb();
    double y=-log(model_pointer->p(ale_pointer));
    //cout <<endl<< "delta=" << delta << "\t tau=" << tau << "\t lambda=" << lambda << "\t ll=" << -y <<endl;
    fval_ = y;
  }
};


int main(int argc, char ** argv)
{
  cout << "ALEml using ALE v"<< ALE_VERSION <<endl;

  if (argc<3) 
    {
      cout << "usage:\n ./ALEml species_tree.newick gene_tree_sample.ale" << endl;
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

  int D=3;
  if (argc>3)
    D=atoi(argv[3]);


  model->set_model_parameter("min_D",D);
  model->set_model_parameter("grid_delta_t",0.005);

  model->construct(Sstring);
  model->set_model_parameter("event_node",0);

  //a set of inital rates
  scalar_type delta=0.01,tau=0.01,lambda=0.1;  
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  //calculate_EGb() must always be called after chaging rates to calculate E-s and G-s 
  //cf. http://arxiv.org/abs/1211.4606
  model->calculate_EGb();

  cout << "Reconciliation model initialised, starting DTL rate optimisation" <<".."<<endl;

    
  //we use the Nelder–Mead method implemented in Bio++ 
  Function* f = new p_fun(model,ale,delta,tau,lambda);
  Optimizer* optimizer = new DownhillSimplexMethod(f);

  optimizer->setProfiler(0);
  optimizer->setMessageHandler(0);
  optimizer->setVerbose(2);
  
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  
  optimizer->init(f->getParameters()); //Here we optimizer all parameters, and start with the default values.
  //FunctionStopCondition stop(optimizer, 1e-1);
  //optimizer->setStopCondition(stop);
  
  optimizer->optimize();
  //optimizer->getParameters().printParameters(cout);
  delta=optimizer->getParameterValue("delta");
  tau=optimizer->getParameterValue("tau");
  lambda=optimizer->getParameterValue("lambda");  
  scalar_type mlll=-optimizer->getFunctionValue();
  cout << endl << "ML rates: " << " delta=" << delta << "; tau=" << tau << "; lambda="<<lambda<<"."<<endl;
  cout << "LL=" << mlll << endl;

  cout << "Calculating ML reconciled gene tree.."<<endl;
 
  //a finer grained reconciliation model for recovering the ML reconciliation
  exODT_model* ml_model=new exODT_model();
  ml_model->set_model_parameter("D",10);
  ml_model->construct(Sstring);
  ml_model->set_model_parameter("delta", delta);
  ml_model->set_model_parameter("tau", tau);
  ml_model->set_model_parameter("lambda", lambda);
  ml_model->calculate_EGb();
  cout << "Reconciliation model initialised" <<".."<<endl;

  //ML backtracking is not compatible with the event node approximation
  ml_model->set_model_parameter("event_node",0);
  //we want events on the leaves as well
  ml_model->set_model_parameter("leaf_events",1);
  //we trace back along the sum to recover the ML reconciled gene tree
  pair<string, scalar_type> res = ml_model->p_MLRec(ale);    
  //and output it..
  string outname=ale_file+".ml_rec"; 
  ofstream fout( outname.c_str() );
  fout <<  "#ALEml using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl<<endl;
  fout << "S:\t"<<ml_model->string_parameter["S_with_ranks"] <<endl;
  fout << endl;
  fout << "Input ale from:\t"<<ale_file<<endl;
  fout << "rate of\t Duplications\tTransfers\tLosses" <<endl;
  fout << "ML \t"<< delta << "\t" << tau << "\t" << lambda << endl;
  fout << endl;

  fout << "reconciled G:\t"<< res.first <<endl;
  fout << endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
  fout <<"Total \t"<< ml_model->MLRec_events["D"] << "\t" << ml_model->MLRec_events["T"] << "\t" << ml_model->MLRec_events["L"]<< "\t" << ml_model->MLRec_events["S"] <<endl;    
  fout << endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tgene copies" <<endl; 
  fout << ml_model->counts_string();
  
  cout << "Results in: " << outname << endl;
  return 1;
}

