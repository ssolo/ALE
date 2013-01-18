#pragma once

#define ALE_VERSION "0.1"

#include <iostream>
#include <Bpp/Numeric/Random.all>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTools.h>
#include <set>
#include <boost/progress.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/algorithm/string.hpp>

//#include <boost/mpi.hpp>

typedef bpp::TreeTemplate<bpp::Node> tree_type;
typedef long double scalar_type;
typedef std::pair<bpp::Node*, bpp::Node*> dedge_type;
//typedef std::pair<long int,std::set < long int > > dip_type;

class approx_posterior
{ 
  
 public:
  //must load
  scalar_type observations;


  //no need to load
  std::string constructor_string;  
  scalar_type alpha;
  scalar_type beta;
  std::map <std::string,int> tree_counts;

  ~approx_posterior()
    {
      leaf_ids.clear();
      id_leaves.clear();
      set_ids.clear(); 
      id_sets.clear();
      Bip_counts.clear();
      for (std::map <long int,std::map< std::set<long int> ,scalar_type> >::iterator it=Dip_counts.begin();it!=Dip_counts.end();it++)    
	{
	  (*it).second.clear();
	}
      Dip_counts.clear();
      Gamma.clear();
      tree_bipstrings.clear();
      bipstring_trees.clear();
      set_sizes.clear();
      for (std::map <int, std::vector <long int > > :: iterator it = size_ordered_bips.begin(); it != size_ordered_bips.end(); it++)
	(*it).second.clear();
      size_ordered_bips.clear();
      Bip_bls.clear();
    }

  approx_posterior();
  approx_posterior(std::string tree);
  void construct(std::string tree_string);

  void save_state(std::string fname);
  void load_state(std::string fname);

  void observation(std::vector<std::string> trees,bool count_topologies=false);
  scalar_type p(std::string tree);

  scalar_type binomial(int n,int m);
  scalar_type trinomial(int n1,int n2,int n3);

  std::pair<std::string,scalar_type> mpp_tree();
  std::string mpp_backtrack(long int g_id, std::map<long int, scalar_type > * qmpp);
  std::string random_tree();
  std::vector<std::string>  all_trees(){return all_trees(Gamma);};//del-loc

  //no need to load
  scalar_type K_Gamma;
  scalar_type N_Gamma;
  std::string name_seperator;
  std::set<int> Gamma;//del-loc
  int Gamma_size;
  std::map <std::string,std::string> tree_bipstrings;//del-loc
  std::map <std::string,std::string> bipstring_trees;//del-loc
  std::map <long int,int> set_sizes;//del-loc
  std::map <int, std::vector <long int > > size_ordered_bips;//del-loc

  //must load
  long int last_leafset_id;
  std::map<std::string,int> leaf_ids;//del-loc
  std::map<int,std::string> id_leaves;//del-loc
  std::map <long int,scalar_type> Bip_counts;//del-loc
  std::map <long int,scalar_type> Bip_bls;//del-loc

  std::map <long int,std::map< std::set<long int>,scalar_type> > Dip_counts;//del-loc
  std::map <std::set <int>,long int>  set_ids;//del-loc  
  std::map< long int, std::set <int> > id_sets;//del-loc
  //nuisance vars
  boost::timer * t;

  //algorithmic
  void decompose(std::string G_string,std::set<int> * bip_ids=NULL );
  std::map <std::set<int>,scalar_type> recompose(std::string G_string);
  void register_leafset(std::string);
  long int set2id(std::set<int> leaf_set);    

  //numeric
  scalar_type Tri(int n2,int n3);
  scalar_type Bi(int n2);
  
  scalar_type p_dip(long int g_id,long int gp_id,long int gpp_id);
  scalar_type p_dip(std::set<int> gamma,std::set<int> gammap,std::set<int> gammapp);
  scalar_type p_bip(long int g_id);
  scalar_type p_bip(std::set<int> gamma);

  //nuisance
  std::string set2name(std::set<int> leaf_set);
  std::string random_split(std::set<int> gamma);
  std::vector<std::string>  all_trees(std::set <int > gamma);//del-loc
 private:
  ;
};

//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
template<typename Set> std::set<Set> powerset(const Set& s, size_t n);
//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
template<typename Set> std::set<Set> powerset(const Set& s);
