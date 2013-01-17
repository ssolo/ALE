
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
    //cout << "#EGb.." << endl;
    double y=-log(model_pointer->p(ale_pointer));
    cout <<endl<< "? d=" << delta << "\t t=" << tau << "\t l=" << lambda << "\t -ll=" << y <<endl;
    fval_ = y;
  }
};


int main(int argc, char ** argv)
{
  cout << "ML ALE v0.1" <<endl;
  string Sstring;
  ifstream file_stream_S (argv[1]);
  getline (file_stream_S,Sstring);

  exODT_model* model=new exODT_model();
  model->set_model_parameter("D",1);
  model->set_model_parameter("DD",10);
  model->construct(Sstring);
  model->set_model_parameter("event_node",1);

  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);

  vector <string> tokens;
  boost::split(tokens,ale_file,boost::is_any_of("/"),boost::token_compress_on);
  ale_file=tokens[tokens.size()-1];


  scalar_type delta=0.01;
  scalar_type tau=0.01;
  scalar_type lambda=0.1;
  
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->calculate_EGb();
  cout << "..1" <<endl;

        
  Function* f = new p_fun(model,ale,delta,tau,lambda);
  Optimizer* optimizer = new DownhillSimplexMethod(f);
  //Optimizer* optimizer = new DownhillSimplexMethod(f);
  //Optimizer* optimizer = new SimpleMultiDimensions(f);

  optimizer->setProfiler(0);
  optimizer->setMessageHandler(0);
  optimizer->setVerbose(0);
  
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  
  optimizer->init(f->getParameters()); //Here we optimizer all parameters, and start with the default values.
  //FunctionStopCondition stop(optimizer, 1e-1);
  //optimizer->setStopCondition(stop);
  cout << "..2" <<endl;
  
  optimizer->optimize();
  //optimizer->getParameters().printParameters(cout);
  delta=optimizer->getParameterValue("delta");
  tau=optimizer->getParameterValue("tau");
  lambda=optimizer->getParameterValue("lambda");  

  exODT_model* max_model=new exODT_model();
  max_model->set_model_parameter("D",10);
  max_model->set_model_parameter("DD",10);
  max_model->construct(Sstring);
  max_model->set_model_parameter("event_node",0);
  max_model->set_model_parameter("delta", delta);
  max_model->set_model_parameter("tau", tau);
  max_model->set_model_parameter("lambda", lambda);

  max_model->set_model_parameter("leaf_events",1);
  max_model->calculate_EGb();

  pair<string, scalar_type> res = max_model->p_MLRec(ale);    

  cout << res.first << endl;    
  string outname=ale_file+".ml_rec"; 
  ofstream fout( outname.c_str() );
  fout << "ML ALE version 0.1 by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl<<endl;
  fout << "S:\t"<<max_model->string_parameter["S_with_ranks"] <<endl;
  fout << endl;
  fout << "Input ale from:\t"<<ale_file<<endl;
  fout << "rate of\t Duplications\tTransfers\tLosses" <<endl;
  fout << "ML \t"<< delta << "\t" << tau << "\t" << lambda << endl;
  fout << endl;

  fout << "reconciled G:\t"<< res.first <<endl;
  fout << endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
  fout <<"Total \t"<< max_model->MLRec_events["D"] << "\t" << max_model->MLRec_events["T"] << "\t" << max_model->MLRec_events["L"]<< "\t" << max_model->MLRec_events["S"] <<endl;    
  fout << endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tgene copies" <<endl; 
  fout << max_model->counts_string();
  
  
  return 1;
  

}

