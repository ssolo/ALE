#include "ALE_util.h"
#include "mpi_tree.h"

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>

using namespace std;
using namespace bpp;
using namespace boost::mpi;

class p_fun:
  public virtual Function,
  public AbstractParametrizable
{
private:
  double fval_;
  mpi_tree* model_pointer;
  int last_branch;
  communicator world;

public:
  p_fun(mpi_tree* model, int last_branch_in,communicator world_in ,double delta_start=0.01,double tau_start=0.01,double lambda_start=0.01) : AbstractParametrizable(""), fval_(0), model_pointer(model) 
  {
    last_branch=last_branch_in;
    world=world_in;
    //We declare parameters here:
    //   IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
    IntervalConstraint* constraint = new IntervalConstraint ( 1e-10, 1e3, true, true );
    for (int e=0;e<last_branch;e++)
	{
	  stringstream deltae;
	  deltae<<"rm_delta"<<"_"<<e;
	  addParameter_( new Parameter(deltae.str(), 1, constraint) ) ;

	  stringstream tau_to_e;
	  tau_to_e<<"rm_tau_to"<<"_"<<e;	  
	  addParameter_( new Parameter(tau_to_e.str(), 1, constraint) ) ;

	  stringstream tau_from_e;
	  tau_from_e<<"rm_tau_from"<<"_"<<e;	  
	  addParameter_( new Parameter(tau_from_e.str(), 1, constraint) ) ;

	  stringstream lambdae;
	  lambdae<<"rm_lambda"<<"_"<<e;
	  addParameter_( new Parameter(lambdae.str(), 1, constraint) ) ;

	  stringstream Oe;
	  Oe<<"rm_O"<<"_"<<e;
	  addParameter_( new Parameter(Oe.str(), 1, constraint) ) ;

	}
    addParameter_( new Parameter("delta", delta_start, constraint) ) ;
    addParameter_( new Parameter("tau", tau_start, constraint) ) ;
    addParameter_( new Parameter("lambda", lambda_start, constraint) ) ;   

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
      
    for (int e=0;e<last_branch;e++)
	{
	  stringstream deltae;
	  deltae<<"rm_delta"<<"_"<<e;
	  double tmp = getParameterValue(deltae.str());
	  model_pointer->model->vector_parameter["rate_multiplier_delta"][e]=tmp;

	  stringstream tau_to_e;
	  tau_to_e<<"rm_tau_to"<<"_"<<e;	  
	  tmp = getParameterValue(tau_to_e.str());	  
	  model_pointer->model->vector_parameter["rate_multiplier_tau_to"][e]=tmp;

	  stringstream tau_from_e;
	  tau_from_e<<"rm_tau_from"<<"_"<<e;
	  tmp = getParameterValue(tau_from_e.str());	  
	  model_pointer->model->vector_parameter["rate_multiplier_tau_from"][e]=tmp;

	  stringstream lambdae;
	  lambdae<<"rm_lambda"<<"_"<<e;
	  tmp = getParameterValue(lambdae.str());
	  model_pointer->model->vector_parameter["rate_multiplier_lambda"][e]=tmp;

	  stringstream Oe;
	  Oe<<"rm_O"<<"_"<<e;
	  tmp = getParameterValue(Oe.str());
	  model_pointer->model->vector_parameter["rate_multiplier_O"][e]=tmp;

	}
    
        double delta = getParameterValue("delta");
        double tau = getParameterValue("tau");
        double lambda = getParameterValue("lambda");

        
        model_pointer->model->set_model_parameter("delta",delta);
        model_pointer->model->set_model_parameter("tau",tau);
        model_pointer->model->set_model_parameter("lambda",lambda);
	

        double y=-(model_pointer->calculate_pun());
        if (world.rank()==0) {
	  cout <<endl<< "delta=" << delta << "\t tau=" << tau << "\t lambda=" << lambda << "\t ll=" << -y <<endl;
	    };
        fval_ = y;
    }
};


int main(int argc, char ** argv)
{

  environment env(argc, argv);
  communicator world;
  int done=1;

  ifstream file_stream_S (argv[1]);
  string Sstring;
  
  getline (file_stream_S,Sstring);
  map<string,scalar_type> parameters;
  if (world.rank()==0) cout << Sstring << endl;
  mpi_tree * infer_tree = new mpi_tree(Sstring,world,parameters,true);
  
  if (world.rank()==0) cout << "..construct.. " << endl;

  infer_tree->load_distributed_ales(argv[2]);
  if (world.rank()==0) cout << "..load.. " << endl;

  
  if (world.rank()==0) cout << infer_tree->model->string_parameter["S_with_ranks"] << endl;
  infer_tree->gather_counts();
  infer_tree->gather_T_to_from();
  broadcast(world,done,0);
  scalar_type delta=0.1;
  scalar_type tau=0.1;
  scalar_type lambda=0.1;
  scalar_type O_R=1,beta=1;
  infer_tree->model->set_model_parameter("O_R",O_R);
  infer_tree->model->set_model_parameter("seq_beta",beta);

  if (atoi(argv[3])==1) infer_tree->model->set_model_parameter("undatedBL",true);
  else infer_tree->model->set_model_parameter("undatedBL",false);
  
  infer_tree->model->calculate_undatedEs();
  
  scalar_type samples=1;//atoi(argv[4]);

  if (argc<7)
    {
  
      Function* f = new p_fun(infer_tree,infer_tree->model->last_branch,world);
      Optimizer* optimizer = new DownhillSimplexMethod(f);
      
      optimizer->setProfiler(0);
      optimizer->setMessageHandler(0);
      optimizer->setVerbose(0);
      if (world.rank()==0)   optimizer->setVerbose(1);
      
      
      optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      cout << "optimizer up to here" <<endl;

      optimizer->init(f->getParameters()); //Here we optimize all parameters, and start with the default values.

      if (world.rank()==0)      cout << "#ML rate optimization.." << endl;
      
      optimizer->optimize();
      broadcast(world,done,0);
      
      for (int e=0;e<infer_tree->model->last_branch;e++)
	{
	  stringstream deltae;
	  deltae<<"rm_delta"<<"_"<<e;
	  double tmp = optimizer->getParameterValue(deltae.str());
	  infer_tree->model->vector_parameter["rate_multiplier_delta"][e]=tmp;
	  
	  stringstream tau_to_e;
	  tau_to_e<<"rm_tau_to"<<"_"<<e;	  
	  tmp = optimizer->getParameterValue(tau_to_e.str());	  
	  infer_tree->model->vector_parameter["rate_multiplier_tau_to"][e]=tmp;
	  
	  stringstream tau_from_e;
	  tau_from_e<<"rm_tau_from"<<"_"<<e;
	  tmp = optimizer->getParameterValue(tau_from_e.str());	  
	  infer_tree->model->vector_parameter["rate_multiplier_tau_from"][e]=tmp;
	  
	  stringstream lambdae;
	  lambdae<<"rm_lambda"<<"_"<<e;
	  tmp = optimizer->getParameterValue(lambdae.str());
	  infer_tree->model->vector_parameter["rate_multiplier_lambda"][e]=tmp;
	  
	  stringstream Oe;
	  Oe<<"rm_O"<<"_"<<e;
	  tmp = optimizer->getParameterValue(Oe.str());
	  infer_tree->model->vector_parameter["rate_multiplier_O"][e]=tmp;
	  
	}
      
      delta = optimizer->getParameterValue("delta");
      tau = optimizer->getParameterValue("tau");
      lambda = optimizer->getParameterValue("lambda");
      
      
      infer_tree->model->set_model_parameter("delta",delta);
      infer_tree->model->set_model_parameter("tau",tau);
      infer_tree->model->set_model_parameter("lambda",lambda);
      if (world.rank()==0 )
	{
	  optimizer->getParameters().printParameters(cout);
	  cout <<endl<< delta << " " << tau << " " << lambda// << " " << sigma
	       << endl;
	}
      
    }
  else
    {
      if (world.rank()==0) cout << "#skipping with: delta=" <<delta <<" lambda="<<lambda<<" tau="<<tau<<endl;
    }
  //optimizer->getParameters().printParameters(cout);
  if (argc>7)
    delta=atof(argv[5]),tau=atof(argv[6]),lambda=atof(argv[7]);samples=atoi(argv[8]);  
  if (world.rank()==0) cout << "#rates : delta=" <<delta <<" lambda="<<lambda<<" tau="<<tau<<endl;


  
  infer_tree->model->set_model_parameter("delta",delta);
  infer_tree->model->set_model_parameter("tau",tau);
  infer_tree->model->set_model_parameter("lambda",lambda);
	

  
  infer_tree->calculate_pun();
  infer_tree->calculate_pun();
  infer_tree->calculate_pun();
  infer_tree->calculate_pun();
  infer_tree->calculate_pun();
  infer_tree->calculate_pun();
  infer_tree->calculate_pun();

  samples=1;
  if (world.rank()==0) cout << "#sampling .." << endl;       
  scalar_type ll_final = infer_tree->calculate_pun(samples);
  if (world.rank()==0) cout << "#sampling done." << endl;       

  infer_tree->gather_counts(samples);
  if (world.rank()==0) cout << "#gather done." << endl;       

  infer_tree->gather_T_to_from(samples);
  if (world.rank()==0) cout << "#gather T_from done." << endl;       

  if (world.rank()==0) cout<< ">tree:\t"<< infer_tree->model->string_parameter["S_with_ranks"] << endl;
  if (world.rank()==0) cout<< ">logl:\t"<< ll_final << endl;
  if (world.rank()==0) cout<< ">Ts:\tfrom\tto"<< endl;
  if (world.rank()==0) infer_tree->print_branch_counts(samples);
  return 0;
  if (world.rank()==0)
    for (map <scalar_type, vector< int > >::iterator it=infer_tree->sort_e.begin();it!=infer_tree->sort_e.end();it++)
      {
	scalar_type Ts=-(*it).first;
	if (Ts>0)
	  for (int i=0;i<(*it).second.size();i++)
	    {
	      int e=infer_tree->sort_e[-Ts][i];
	      int f=infer_tree->sort_f[-Ts][i];
	      if (e<infer_tree->model->last_leaf)
		cout << "\t" << infer_tree->model->node_name[infer_tree->model->id_nodes[e]];
	      else
		cout << "\t" << e;
	      if (f<infer_tree->model->last_leaf)
		cout << "\t" << infer_tree->model->node_name[infer_tree->model->id_nodes[f]];	      
	      else
		cout << "\t" << f;
	      cout << "\t" << Ts << endl; //" " << new_S << endl;		
	    }
      }


}
