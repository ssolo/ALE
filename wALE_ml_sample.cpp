
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
    //   IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
    IntervalConstraint* constraint = new IntervalConstraint ( 1e-6, 10-1e-6, true, true );
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
    double delta = getParameterValue("delta");
    double tau = getParameterValue("tau");
    double lambda = getParameterValue("lambda");
    model_pointer->set_model_parameter("delta",delta);
    model_pointer->set_model_parameter("tau",tau);
    model_pointer->set_model_parameter("lambda",lambda);
    model_pointer->calculate_EGb();
    double y=-log(model_pointer->p(ale_pointer));
    cout <<endl<< "delta=" << delta << "\t tau=" << tau << "\t lambda=" << lambda << "\t ll=" << -y <<endl;
    fval_ = y;
  }
};


int main(int argc, char ** argv)
{
  cout << "ALEml using ALE v"<< ALE_VERSION <<endl;

  if (argc<3) 
    {
      cout << "usage:\n ./ALEml species_tree.newick gene_tree_sample.ale [gene_name_seperator]" << endl;
      return 1;
    }

  //we need a dated species tree in newick format
  string Sstring;
  ifstream file_stream_S (argv[1]);
  string S_tree_name=argv[1];
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
  cout << "o" << endl;

  int D=4;
  model->set_model_parameter("gene_name_separators", ".");
  model->set_model_parameter("BOOT_STRAP_LABLES","yes");

  model->set_model_parameter("min_D",D);
  model->set_model_parameter("grid_delta_t",0.05);

  model->construct(Sstring);
  model->set_model_parameter("event_node",0);
  model->set_model_parameter("leaf_events",1);
  model->set_model_parameter("N",1);

  //a set of inital rates
  scalar_type delta=0.01,tau=0.01,lambda=0.1;  
  if (argc>6)
    delta=atof(argv[4]),tau=atof(argv[5]),lambda=atof(argv[6]);  
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("sigma_hat", 1);

  //calculate_EGb() must always be called after changing rates to calculate E-s and G-s
  //cf. http://arxiv.org/abs/1211.4606
  model->calculate_EGb();

  cout << "Reconciliation model initialised, starting DTL rate optimisation" <<".."<<endl;
  
  if (true)
    {
      
      //we use the Nelder–Mead method implemented in Bio++
      Function* f = new p_fun(model,ale,delta,tau,lambda);
      Optimizer* optimizer = new DownhillSimplexMethod(f);
      
      optimizer->setProfiler(0);
      optimizer->setMessageHandler(0);
      optimizer->setVerbose(2);
      
      optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      optimizer->init(f->getParameters()); //Here we optimize all parameters, and start with the default values.
  
      //   FunctionStopCondition stop(optimizer, 1);//1e-1);
      // optimizer->setStopCondition(stop);
      //TEMP
      //optimizer->setMaximumNumberOfEvaluations( 10 );
      
      optimizer->optimize();
      
      //optimizer->getParameters().printParameters(cout);
      delta=optimizer->getParameterValue("delta");
      tau=optimizer->getParameterValue("tau");
      lambda=optimizer->getParameterValue("lambda");

      scalar_type mlll=-optimizer->getFunctionValue();
      cout << endl << "ML rates: " << " delta=" << delta << "; tau=" << tau << "; lambda="<<lambda<<"."<<endl;
      cout << "LL=" << mlll << endl;

      cout << "Calculating ML reconciled gene tree.."<<endl;
 

      pair<string, scalar_type> res = model->p_MLRec(ale);    
      //and output it..
      string outname=S_tree_name+".ml_rec"; 

      ofstream fout( outname.c_str() );
      fout <<  "#ALEml using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl<<endl;
      fout << "S:\t"<<model->string_parameter["S_with_ranks"] <<endl;
      fout << endl;
      fout << "Input ale from:\t"<<ale_file<<endl;
      fout << "rate of\t Duplications\tTransfers\tLosses" <<endl;
      fout << "ML \t"<< delta << "\t" << tau << "\t" << lambda << endl;
      fout << endl;

      fout << "reconciled G:\t"<< res.first <<endl;
      fout << endl;
      fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
      fout <<"Total \t"<< model->MLRec_events["D"] << "\t" << model->MLRec_events["T"] << "\t" << model->MLRec_events["L"]<< "\t" << model->MLRec_events["S"] <<endl;    
      fout << endl;
      fout << "# of\t Duplications\tTransfers\tLosses\tgene copies" <<endl; 
      fout << model->counts_string();
  
      cout << "Results in: " << outname << endl;

    }

  cout << "strating sampling.. " << endl;
  string sample_name=S_tree_name+".ml_samples"; 
  ofstream sample_out( sample_name.c_str() );

  string hist_name=S_tree_name+".hist"; 

  map <string,map <int,map<int,int> > > event_histograms;
  vector<string> event_types;
  event_types.push_back("Os");
  event_types.push_back("Ds");
  event_types.push_back("Ts");
  event_types.push_back("Ls");
  event_types.push_back("Tfroms");
  event_types.push_back("copies");
  scalar_type samples_in_histogram=0;

  for (vector<string>::iterator et=event_types.begin();et!=event_types.end();++et )
    for (int branch=0;branch<model->last_branch;branch++)
      {
	for (int i=0;i<100;++i)
	  event_histograms[(*et)][branch][i]=0;
	//model->branch_counts[(*et)][branch];      
      }	      

  samples_in_histogram++;
  for (vector<string>::iterator et=event_types.begin();et!=event_types.end();++et )		
    for (int branch=0;branch<model->last_branch;branch++)
      {
	event_histograms[(*et)][branch][ model->branch_counts[(*et)][branch] ]++;		  
	model->branch_counts[(*et)][branch]=0;
      }

  cout << model->p(ale) << endl;
  
  for (int i =0; i < 10000; i++ )
    {
      string sample_tree=model->sample(false);
      //cout << sample_tree << endl;
      samples_in_histogram++;
      for (vector<string>::iterator et=event_types.begin();et!=event_types.end();++et )		
	for (int branch=0;branch<model->last_branch;branch++)
	  {
	    event_histograms[(*et)][branch][ model->branch_counts[(*et)][branch] ]++;		  
	    model->branch_counts[(*et)][branch]=0;
	  }
      sample_out << sample_tree << endl;
    }
  cout << "sampling done." << endl;
  for (vector<string>::iterator et=event_types.begin();et!=event_types.end();++et )
    {
      ofstream h_out( S_tree_name+"_"+((*et)+".h").c_str() );
      h_out << "# Fraction of samples with number of " << (*et) << " in " <<samples_in_histogram << " samples for each branch of S." << endl;
      h_out << "#id\ttb\tte\trnk";
      for (int i=0;i<10;++i)
	h_out << "\t" << i;
      h_out<<endl;

      for (int branch=0;branch<model->last_branch;branch++)
	{

	  stringstream named_branch;
	  if (branch==model->alpha)
	    named_branch<<-1;
	  else if (model->id_ranks[branch]==0)
	    named_branch<<model->extant_species[branch];
	  else
	    {
	      if (model->rank2label[model->id_ranks[branch]]!=-1)
		named_branch<<model->rank2label[model->id_ranks[branch]];
	      else
		named_branch<<"ROOT";
	    }

	  h_out<<branch;
	  h_out<<"\t"<<model->t_begin[branch];
	  h_out<<"\t"<<model->t_end[branch];
	  h_out<<"\t"<<named_branch.str();
	  for (int i=0;i<10;++i)
	    h_out << "\t" << event_histograms[(*et)][branch][i]/samples_in_histogram;
	  h_out<<endl;
	  //model->branch_counts[(*et)][branch];      
	}	      
    }

  return 0;
}

