#include "ALE_util.h"
#include "mpi_tree.h"

#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>

using namespace std;
using namespace bpp;
using namespace boost::mpi;

class p_fun : public virtual Function, public AbstractParametrizable {
private:
  double fval_;
  mpi_tree *model_pointer;

public:
  p_fun(mpi_tree *model, double delta_start = 0.2, double tau_start = 0.2,
        double lambda_start = 0.5 //,double sigma_start=2
        )
      : AbstractParametrizable(""), fval_(0), model_pointer(model) {
    // We declare parameters here:
    //   IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
    IntervalConstraint *constraint =
        new IntervalConstraint(1e-6, 10 - 1e-6, true, true);
    addParameter_(new Parameter("delta", delta_start, constraint));
    addParameter_(new Parameter("tau", tau_start, constraint));
    addParameter_(new Parameter("lambda", lambda_start, constraint));
    // addParameter_( new Parameter("sigma", sigma_start, constraint) ) ;
  }

  p_fun *clone() const { return new p_fun(*this); }

public:
  void setParameters(const ParameterList &pl) throw(ParameterNotFoundException,
                                                    ConstraintException,
                                                    Exception) {
    matchParametersValues(pl);
  }
  double getValue() const throw(Exception) { return fval_; }
  void fireParameterChanged(const ParameterList &pl) {
    double delta = getParameterValue("delta");
    double tau = getParameterValue("tau");
    double lambda = getParameterValue("lambda");

    // double sigma = getParameterValue("sigma");

    model_pointer->model->set_model_parameter("delta", delta);
    model_pointer->model->set_model_parameter("tau", tau);
    model_pointer->model->set_model_parameter("lambda", lambda);

    // model_pointer->model->set_model_parameter("Delta_bar",sigma*1e6);
    // model_pointer->model->set_model_parameter("Lambda_bar",sigma*1e6);

    // model_pointer->calculate_EGb();
    double y = -(model_pointer->calculate_p());
    // if (world.rank()==0) cout <<endl<< "delta=" << delta << "\t tau=" << tau
    // << "\t lambda=" << lambda << "\t ll=" << -y <<endl;
    fval_ = y;
  }
};

int main(int argc, char **argv) {
  environment env(argc, argv);
  communicator world;
  vector<string> ale_names;
  if (world.rank() == 0) {
    ifstream file_stream(argv[2]);
    while (!file_stream.eof()) {
      string fname;
      getline(file_stream, fname);
      boost::trim(fname);
      if (fname.find(".ale") != fname.npos)
        ale_names.push_back(fname);
    }
  }
  int done = 1;
  if (world.rank() == 0)
    cout << "#list of " << ale_names.size() << " ales.." << endl;
  broadcast(world, done, 0);

  ifstream file_stream_S(argv[1]);
  string Sstring;
  getline(file_stream_S, Sstring);
  mpi_tree *infer_tree = new mpi_tree(Sstring, world);

  infer_tree->distribute_ales(ale_names);
  scalar_type ll = infer_tree->calculate_p();
  if (world.rank() == 0)
    cout << "LL = " << ll << endl;
  broadcast(world, done, 0);

  Function *f = new p_fun(infer_tree);
  Optimizer *optimizer = new DownhillSimplexMethod(f);

  optimizer->setProfiler(0);
  optimizer->setMessageHandler(0);
  optimizer->setVerbose(0);
  if (world.rank() == 0)
    optimizer->setVerbose(0);

  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  optimizer->init(f->getParameters()); // Here we optimize all parameters, and
                                       // start with the default values.

  //   FunctionStopCondition stop(optimizer, 1);//1e-1);
  // optimizer->setStopCondition(stop);
  // TEMP
  // optimizer->setMaximumNumberOfEvaluations( 10 );

  optimizer->optimize();

  // optimizer->getParameters().printParameters(cout);
  broadcast(world, done, 0);

  if (world.rank() == 0) {
    optimizer->getParameters().printParameters(cout);
    scalar_type delta = optimizer->getParameterValue("delta");
    scalar_type tau = optimizer->getParameterValue("tau");
    scalar_type lambda = optimizer->getParameterValue("lambda");
    // scalar_type sigma=optimizer->getParameterValue("sigma");

    cout << endl
         << delta << " " << tau << " " << lambda // << " " << sigma
         << endl;
  }
}
