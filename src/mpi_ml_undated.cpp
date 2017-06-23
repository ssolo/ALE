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
  communicator world;
public:
  p_fun(mpi_tree* model, communicator world_in , double delta_start=0.05,double tau_start=0.,double lambda_start=0.2//,double sigma_start=2
  ) : AbstractParametrizable(""), fval_(0), model_pointer(model)
  {
    world=world_in;

    //We declare parameters here:
    //   IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
    IntervalConstraint* constraint = new IntervalConstraint ( 0, 4, true, true );
    addParameter_( new Parameter("delta", delta_start, constraint) ) ;
    addParameter_( new Parameter("tau", tau_start, constraint) ) ;
    addParameter_( new Parameter("lambda", lambda_start, constraint) ) ;
    //addParameter_( new Parameter("sigma", sigma_start, constraint) ) ;

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

    //double sigma = getParameterValue("sigma");

    model_pointer->model->set_model_parameter("delta",delta);
    model_pointer->model->set_model_parameter("tau",tau);
    model_pointer->model->set_model_parameter("lambda",lambda);

    //model_pointer->model->set_model_parameter("Delta_bar",sigma*1e6);
    //model_pointer->model->set_model_parameter("Lambda_bar",sigma*1e6);

    //model_pointer->calculate_EGb();
    double y=-(model_pointer->calculate_pun());
    //model_pointer->gather_counts();
    if (world.rank()==0) {
      cout <<endl<< "delta=" << delta << "\t tau=" << tau << "\t lambda=" << lambda << "\t ll=" << -y <<endl;
      //model_pointer->print_branch_counts();
    };
    fval_ = y;
  }
};


int main(int argc, char ** argv)
{

  environment env(argc, argv);
  communicator world;
  int done=1;

  if (argc<2)
  {
    cout << "usage:\n ./mpi_ml_undated species_tree_file file_containing_ales O_R=OriginationAtRoot delta=DuplicationRate tau=TransferRate lambda=LossRate beta=weight_of_sequence_evidence outputFiles=n" << endl;
    cout << "\n\tImportant: If duplication, transfer or loss rates are specified, they are fixed." << endl;
    return 1;
  }


  if (!fexists(argv[1])) {
    cout << "Error, file "<<argv[1] << " does not seem accessible." << endl;
    exit(1);
  }
  ifstream file_stream_S (argv[1]);
  string Sstring;

  getline (file_stream_S,Sstring);
  map<string,scalar_type> parameters;
  if (world.rank()==0) cout << Sstring << endl;
  mpi_tree * infer_tree = new mpi_tree(Sstring,world,parameters,true);
  if (world.rank()==0) cout << "Constructed the species tree. " << endl;

  infer_tree->load_distributed_ales(argv[2]);
  if (world.rank()==0) cout << "Ales loaded. " << endl;

  //scalar_type ll = infer_tree->calculate_pun(3);

  //scalar_type ll = infer_tree->calculate_pun(10,1);

  if (world.rank()==0) cout << infer_tree->model->string_parameter["S_with_ranks"] << endl;
  //infer_tree->gather_counts();
  //infer_tree->gather_T_to_from();
  broadcast(world,done,0);
  scalar_type delta=0.01;
  scalar_type tau=0.01;
  scalar_type lambda=0.1;

  // Getting the other options
  scalar_type samples=100;
  scalar_type O_R=1,beta=1;
  string fractionMissingFile = "";
  bool outputyn = false;
  vector<string> paramsToIgnore;
  for (int i=3;i<argc;i++)
  {
    string next_field=argv[i];
    vector <string> tokens;
    boost::split(tokens,next_field,boost::is_any_of("="),boost::token_compress_on);
    if (tokens[0]=="sample")
    samples=atoi(tokens[1].c_str());
    else if (tokens[0]=="delta")
    {
      delta=atof(tokens[1].c_str());
      cout << "\n\tDelta fixed to " << delta << endl;
      paramsToIgnore.push_back("delta");
    }
    else if (tokens[0]=="tau")
    {
      tau=atof(tokens[1].c_str());
      cout << "\n\tTau fixed to " << tau << endl;
      paramsToIgnore.push_back("tau");
    }
    else if (tokens[0]=="lambda")
    {
      lambda=atof(tokens[1].c_str());
      cout << "Lambda fixed to " << lambda << endl;
      paramsToIgnore.push_back("lambda");
    }
    else if (tokens[0]=="O_R")
    {
      O_R=atof(tokens[1].c_str());
      cout << "\n\tO_R set to " << O_R << endl;
    }
    else if (tokens[0]=="beta")
    {
      beta=atof(tokens[1].c_str());
      cout << "\n\tBeta set to " << beta << endl;
    }
    else if (tokens[0]=="fraction_missing")
    {
      fractionMissingFile=tokens[1];
      cout << "\n\tFile containing fractions of missing genes set to " << fractionMissingFile << endl;
    }
    else if (tokens[0]=="outputFiles")
    {
      if (tokens[1] == "y" || tokens[1] == "yes" || tokens[1] == "Y" || tokens[1] == "YES") {
        outputyn = true;
      }
    }
  }


  Function* f = new p_fun(infer_tree,world );
  Optimizer* optimizer = new DownhillSimplexMethod(f);

  optimizer->setProfiler(0);
  optimizer->setMessageHandler(0);
  optimizer->setVerbose(0);
  if (world.rank()==0)   optimizer->setVerbose(1);


  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  ParameterList parametersRep = f->getParameters();
  //Removing the parameters that were fixed in the input
  for (size_t i=0;i<paramsToIgnore.size();i++) {
    parametersRep.deleteParameter(paramsToIgnore[i]);
  }
  optimizer->init(parametersRep); //Here we optimize all parameters, and start with the default values.

  if ( parametersRep.size() > 0)  {
    if (world.rank()==0) cout << "#ML rate optimization.." << endl;
    optimizer->optimize();
    broadcast(world,done,0);
    delta=optimizer->getParameterValue("delta");
    tau=optimizer->getParameterValue("tau");
    lambda=optimizer->getParameterValue("lambda");
    if (world.rank()==0 )
    {
      optimizer->getParameters().printParameters(cout);
      cout <<endl<< delta << " " << tau << " " << lambda// << " " << sigma
      << endl;
    }
  }
  else {
    if (world.rank()==0) cout << "#skipping rate optimization" << endl;
    if (world.rank()==0) cout << "#rates : delta=" <<delta <<" lambda="<<lambda<<" tau="<<tau<<endl;
  }

  //optimizer->getParameters().printParameters(cout);
  //  if (argc>5)
  //  delta=atof(argv[3]),tau=atof(argv[4]),lambda=atof(argv[5]);
  if (outputyn) {
    infer_tree->model->set_model_parameter("delta",delta);
    infer_tree->model->set_model_parameter("tau",tau);
    infer_tree->model->set_model_parameter("lambda",lambda);
    infer_tree->calculate_pun(0);
    scalar_type samples=1;
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

    string Tname=Sstring+".Ts";
    ofstream fout( Tname.c_str() );

    if (world.rank()==0)
    {
      for (map <scalar_type, vector< int > >::iterator it=infer_tree->sort_e.begin();it!=infer_tree->sort_e.end();it++)
      {
        scalar_type Ts=-(*it).first;
        if (Ts>0)
        for (int i=0;i<(*it).second.size();i++)
        {
          int e=infer_tree->sort_e[-Ts][i];
          int f=infer_tree->sort_f[-Ts][i];
          if (e<infer_tree->model->last_leaf)
          fout << "\t" << infer_tree->model->node_name[infer_tree->model->id_nodes[e]];
          else
          fout << "\t" << e;
          if (f<infer_tree->model->last_leaf)
          fout << "\t" << infer_tree->model->node_name[infer_tree->model->id_nodes[f]];
          else
          fout << "\t" << f;
          fout << "\t" << Ts << endl; //" " << new_S << endl;
        }
      }
    }
    fout.close();
  }
}
