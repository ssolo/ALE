#include "exODT.h"
#include <boost/mpi.hpp>

class mpi_tree
{
 public:
  exODT_model * model;//del-loc
  boost::mpi::communicator world;
  int server,rank,size;
  std::vector<approx_posterior*> ale_pointers;//del-loc
  std::vector<std::string> client_fnames;
  scalar_type N_ales;

  std::vector<std::pair<std::string,scalar_type> > MLRec_res;//del-loc
  std::vector<std::string> Ttokens;//del-loc

  std::map<std::string,scalar_type> scalar_parameter;
  std::map<std::string,std::string> string_parameter;

  void set_parameter(std::string name, scalar_type value)
  {
    scalar_parameter[name]=value;    
  };
  void set_parameter(std::string name, std::string value)
  {
    string_parameter[name]=value;
  };


  scalar_type delta_avg,tau_avg,lambda_avg;
  scalar_type delta_norm,tau_norm,lambda_norm;
  std::vector<scalar_type> delta_branch_avg,tau_branch_avg,lambda_branch_avg;//del-loc
  std::vector<scalar_type> delta_branch_norm,tau_branch_norm,lambda_branch_norm;//del-loc
  std::string S_string;
  mpi_tree(std::string Sstring,const boost::mpi::communicator mpi_world,std::map<std::string,scalar_type> set_parameters=std::map<std::string,scalar_type>(),bool undated=false)
    {
      if (undated)
	{
	  S_string=Sstring;
	  set_parameter("min_delta",1e-6);
	  set_parameter("min_tau",1e-6);
	  set_parameter("min_lambda",1e-6);

	  set_parameter("inital_delta",0.01);
	  set_parameter("inital_tau",0.01);
	  set_parameter("inital_lambda",0.02);      
	  model=new exODT_model();
	  model->construct_undated(Sstring);//del-loc

	  model->set_model_parameter("delta",scalar_parameter["inital_delta"]);
	  model->set_model_parameter("tau",scalar_parameter["inital_tau"]);
	  model->set_model_parameter("lambda",scalar_parameter["inital_lambda"]);
	  for (std::map<std::string,scalar_type>::iterator it=set_parameters.begin();it!=set_parameters.end();it++)
	    model->set_model_parameter((*it).first,(*it).second);

	}
      else
	{
	  set_parameter("use_mpp_trees",0);
	  set_parameter("min_delta",1e-6);
	  set_parameter("min_tau",1e-6);
	  set_parameter("min_lambda",1e-6);
	  
	  set_parameter("inital_delta",0.01);
	  set_parameter("inital_tau",0.01);
	  set_parameter("inital_lambda",0.02);      
	  
	  model=new exODT_model();      
	  model->set_model_parameter("min_D",3);
	  model->set_model_parameter("grid_delta_t",0.005);
	  model->set_model_parameter("event_node",0);
	  model->set_model_parameter("DD",10);
	  for (std::map<std::string,scalar_type>::iterator it=set_parameters.begin();it!=set_parameters.end();it++)
	    model->set_model_parameter((*it).first,(*it).second);
	  model->construct(Sstring);//del-loc
	  scalar_type N=1e6;
	  model->set_model_parameter("N",1e6);//we can almost scale out N assuming height from coalescent..
	  model->set_model_parameter("Delta_bar",N);
	  model->set_model_parameter("Lambda_bar",N);
	  model->set_model_parameter("delta",scalar_parameter["inital_delta"]);
	  model->set_model_parameter("tau",scalar_parameter["inital_tau"]);
	  model->set_model_parameter("lambda",scalar_parameter["inital_lambda"]);
	  
	  model->calculate_EGb();//with default parameters
	}
      world = mpi_world;      
      server=0;
      rank = world.rank();
      size = world.size();


    };
  ~mpi_tree()
    {
      for (std::vector<approx_posterior*>::iterator it=ale_pointers.begin();it!=ale_pointers.end();it++)
	delete (*it);
      ale_pointers.clear();     
      MLRec_res.clear();
      client_fnames.clear();
      delta_branch_avg.clear(),tau_branch_avg.clear(),lambda_branch_avg.clear();//del-loc
      delta_branch_norm.clear(),tau_branch_norm.clear(),lambda_branch_norm.clear();//del-loc

      delete model;
    };

  //implimented in mpi_tree.cpp
  void distribute_ales(std::vector<std::string>,bool list_of_trees=false);
  void load_distributed_ales(std::string fname);

  std::vector < std::vector<scalar_type> > gathered_T_to_from;
  void gather_T_to_from();

  void gather_counts();
  void clear_counts();
  std::string branch_counts_string();
  void show_branch_counts();
  void print_branch_counts();

  scalar_type calculate_MLRecs(bool estimate=false,bool branchwise=false);
  scalar_type calculate_p();
  scalar_type calculate_pun();
  scalar_type calculate_punt(std::string S);
  scalar_type calculate_pun(int n,bool bw=false);
  std::map <scalar_type, std::vector< int > >sort_e;
  std::map <scalar_type, std::vector< int > >sort_f;

  void estimate_rates();
  void estimate_rates_bw();
  
  //implimented in rae_estimate.cpp
  std::vector<scalar_type> dtl_estimate(int branch,scalar_type N_ales_norm);
  //scalar_type estimate_rates(std::string mode="uniform");
  void show_rates();
 private:
  ;
};
