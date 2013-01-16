#include "ALE.h"

struct step 
{
  int e;
  int ep;
  int epp;
  scalar_type t;
  int rank;
  long int g_id;
  long int gp_id;
  long int gpp_id;
  std::string event;
};

class exODT_model
{
 public:
  
  std::map <std::string,scalar_type> scalar_parameter;//del_loc
  std::map <std::string,std::vector <scalar_type> > vector_parameter;//del_loc
  std::map <std::string,std::string> string_parameter;//del_loc

  int signal;
  std::string signal_string;

  int alpha;
  int last_branch;
  int last_rank;
  approx_posterior * ale_pointer;

  std::map<int,int> father;//del-loc
  std::map<int,std::vector<int> > daughters;//del-loc
  std::map<int,std::string> extant_species;//del-loc
  std::map<scalar_type,int> branch_ts;//del-loc
  std::map<int,int>rank_ids;//del-loc
  std::map<int,int>id_ranks;//del-loc

  tree_type * S;
  bpp::Node * S_root;
  std::map<bpp::Node *,int>node_ids;
  std::map<int,bpp::Node *>id_nodes;

  std::map<int,scalar_type > t_begin;//del-loc
  std::map<int,scalar_type > t_end;//del-loc
  
  std::map<int,std::vector<int > > time_slices;//del-loc
  std::map<int,std::vector<int > > branch_slices;//del-loc
  std::map<int,std::vector<scalar_type > > time_slice_times;//del-loc
  std::map<int,scalar_type > time_slice_begins;//del-loc
  
  std::map<int,std::map <scalar_type,scalar_type> > Ee;//del-loc
  std::map<int,std::map <scalar_type,scalar_type> > Ge;//del-loc
  std::map<long int, std::map< scalar_type, std::map<int, scalar_type> > > q;//del-loc
  std::map<long int, std::map< scalar_type, std::map<int, step> > > q_step;//del-loc
  std::map <long int,std::string> gid_sps;//del-loc  
 
  std::map <std::string,scalar_type> MLRec_events;//del-loc  
  std::map<std::string, std::vector<scalar_type> > branch_counts;//del-loc
  std::vector<std::string> Ttokens;//del-loc
  //implimented in exODT.cpp
  void construct(std::string Sstring,scalar_type N=1e6);
  exODT_model();
  ~exODT_model()
    {  
      rank_ids.clear();
      id_ranks.clear();
      father.clear();
      for (std::map<int,std::vector<int> >::iterator it=daughters.begin();it!=daughters.end();it++)
	(*it).second.clear();
      daughters.clear();
      extant_species.clear();
      branch_ts.clear();
      rank_ids.clear();
      id_ranks.clear();
      t_begin.clear();
      t_end.clear();
      for (std::map<int,std::vector<int> >::iterator it=time_slices.begin();it!=time_slices.end();it++)
	(*it).second.clear();
      time_slices.clear();
      for (std::map<int,std::vector<int> >::iterator it=branch_slices.begin();it!=branch_slices.end();it++)
	(*it).second.clear();
      for (std::map<int,std::vector<scalar_type > >::iterator it=time_slice_times.begin();it!=time_slice_times.end();it++)
	(*it).second.clear();
      time_slice_times.clear();
      time_slice_begins.clear();
      scalar_parameter.clear();
      for (std::map <std::string,std::vector <scalar_type> >::iterator it=vector_parameter.begin();it!=vector_parameter.end();it++)//del_loc
	(*it).second.clear();
      vector_parameter.clear();
      string_parameter.clear();
      node_ids.clear();
      id_nodes.clear();
      delete S;
      for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ee.begin();it!=Ee.end();it++)//del_loc
	(*it).second.clear();
      Ee.clear();
      for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ge.begin();it!=Ge.end();it++)//del_loc
	(*it).second.clear();
      Ge.clear();
      Ee.clear();
      for (std::map<long int, std::map< scalar_type, std::map<int, scalar_type> > >::iterator it=q.begin();it!=q.end();it++)
	{
	  for ( std::map< scalar_type, std::map<int, scalar_type> >::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	    (*jt).second.clear();
	  (*it).second.clear();
	}      
      q.clear();
      for (std::map<long int, std::map< scalar_type, std::map<int, step> > >::iterator it=q_step.begin();it!=q_step.end();it++)
	{
	  for ( std::map< scalar_type, std::map<int, step> >::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	    (*jt).second.clear();
	  (*it).second.clear();
	}      
      q_step.clear();
      gid_sps.clear();
      MLRec_events.clear();
      for (std::map<std::string, std::vector<scalar_type> >::iterator it=branch_counts.begin();it!=branch_counts.end();it++)//del_loc
	(*it).second.clear();
      branch_counts.clear();
      Ttokens.clear();
    };

  void set_model_parameter(std::string name,std::string value);
  void set_model_parameter(std::string,scalar_type);
  void set_model_parameter(std::string,std::vector<scalar_type>);

  //implimented in model.cpp
  scalar_type p(approx_posterior *ale);
  void calculate_EG();
  void calculate_EGb();

  //implimented in traceback.cpp
  std::pair<std::string,scalar_type> p_MLRec(approx_posterior *ale,bool lowmem=true);
  std::pair<std::string,scalar_type> traceback();
  std::string traceback(long int g_id,scalar_type t,scalar_type rank,int e,scalar_type branch_length,std::string branch_events,std::string transfer_token="");
  void register_O(int e);
  void register_D(int e);
  void register_Tto(int e);
  void register_Tfrom(int e);
  void register_L(int e);
  void register_S(int e);
  void register_leaf(int e);
  void register_Ttoken(std::string token);
  //implimented in traceback_lowmem.cpp - not done
  std::pair<std::string,scalar_type> p_MLRec_lowmem(approx_posterior *ale);
  std::string traceback_lowmem(long int g_id,scalar_type t,scalar_type rank,int e,scalar_type branch_length,std::string branch_events,std::string transfer_token="");
  //implimented in sample.cpp
  std::string sample(bool max_rec=false);
  std::string sample(bool S_node,long int g_id,int t_i,scalar_type rank,int e,scalar_type branch_length,std::string branch_events, std::string transfer_toke="",bool max_rec=false);


  void show_counts(std::string name);
  std::string counts_string();

  void show_rates(std::string name);


 private:
  ;
  
};
