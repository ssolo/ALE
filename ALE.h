//all code by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
#pragma once

#define ALE_VERSION "0.2"

#include <iostream>
#include <Bpp/Numeric/Random.all>
#include <Bpp/Numeric/NumConstants.h>
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


/****************************************************************************
 // approx_posterior class.
 // This class contains the description of a posterior on phylogenetic trees
 // by its clade probabilities and conditional clade probabilities.
 //Lexicon:
 // leaf set = bipartition = clade
 *****************************************************************************/

class approx_posterior
{ 
  
 public:
  //must load
  scalar_type observations;                    //Number of trees observed to build the approx_posterior object.


  //no need to load
  std::string constructor_string;              //string representing the tree in Newick format
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

  approx_posterior();                                                               //Does nothing. Formal constructor must be followed by load_state.
  approx_posterior(std::string tree);                                               //Constructs a basic instance by calling construct.
  void construct(std::string tree_string);                                          //Constructs a basic instance.

  void save_state(std::string fname);                                               //Writes the object to a file.
  void load_state(std::string fname);

  void observation(std::vector<std::string> trees,bool count_topologies=false);     //Given a vector of trees, fills an approx_posterior object by recursively calling decompose, and then doing some more counts.
  scalar_type p(std::string tree);                                                  //Computes the probability of a string tree. Calls recompose on the tree, and then uses the map returned by recompose to compute the probability of the whole tree.
  scalar_type nbipp(std::string tree);

  scalar_type binomial(int n,int m);                                                //Computes the binomial coefficient.
  scalar_type trinomial(int n1,int n2,int n3);                                      //Computes the multinomial coefficient for 3 elements.

  std::pair<std::string,scalar_type> mpp_tree();                                    //Returns the maximum a posteriori tree that can be amalgamated from the approx_prior object. Uses a double-recursive traversal of all bipartitions.
  std::string mpp_backtrack(long int g_id, std::map<long int, scalar_type > * qmpp);//Recursive function that, given a bipartition id and a map associating bipartition ids to their maximum a posteriori value, builds the maximum a posteriori tree, complete with (average) branch lengths.
  std::string random_tree();                                                        //Function that returns a random tree with unit branch lengths. Calls random_split.
	std::vector<std::string>  all_trees(){return all_trees(Gamma);};                //del-loc. 


  //no need to load
	std::set<int> Gamma;                                                             //del-loc. Set containing all leaf ids. ~Clade of all leaves in the tree.
	int Gamma_size;                                                                  //Number of leaves.
	scalar_type K_Gamma;                                                             //number of bipartitions of Gamma.
  scalar_type N_Gamma;                                                               //number of unrooted trees on Gamma_size leaves
  std::string name_seperator;                                                        //Character used between leaf names when writing the content of a leaf set.
  std::map <std::string,std::string> tree_bipstrings;                                //del-loc. Map between tree string and string containing all bipartitions in the tree.
  std::map <std::string,std::string> bipstring_trees;                                //del-loc. Dual from above. Map between string containing all bipartitions in the tree and tree string.
  std::map <long int,int> set_sizes;                                                 //del-loc. Map between a bipartition id and the sizes of the corresponding leaf set.
  std::map <int, std::vector <long int > > size_ordered_bips;                        //del-loc. Map between bipartition size, and the ids of all bipartitions of this size.

  //must load
  long int last_leafset_id;                                                          //Total number of sets of leaves (=bipartitions) observed in the posterior.
  std::map<std::string,int> leaf_ids;                                                //del-loc. Map between species name and leaf id. Leaf ids go from 1 to Gamma_size.
  std::map<int,std::string> id_leaves;                                               //del-loc. Map between leaf id and species name. Dual from above.
  std::map <long int,scalar_type> Bip_counts;                                        //del-loc. For each bipartition, gives the number of times it was observed.
  std::map <long int,scalar_type> Bip_bls;                                           //del-loc. Sum of the branch lengths associated to the bipartitions.

  std::map <long int,std::map< std::set<long int>,scalar_type> > Dip_counts;         //del-loc. Contains the frequency of triplets: mother clade and its two daughter clades. Map between the bipartition id of the mother clade and another map containing a set of bipartition ids (couldn't it be just a pair, in the case of bifurcating trees?) and the frequency of the associated triplet of bipartitions. 
  std::map <std::set <int>,long int>  set_ids;                                       //del-loc. Map between a set of leaf ids and the corresponding bipartition index. 
  std::map< long int, std::set <int> > id_sets;                                      //del-loc. Dual from above. Map between a bipartition index and the corresponding leaf ids.
  //nuisance vars
  boost::timer * t;

  //algorithmic
  void decompose(std::string G_string,std::set<int> * bip_ids=NULL );                //Parses a tree in string format and updates the approx_prior object accordingly (notably updates the Bip_bls, Bip_counts, Dip_counts, and set_ids + id_sets through set2id)
  std::map <std::set<int>,scalar_type> recompose(std::string G_string);              //For a given input tree string, returns a map between all sets of leaves contained in the tree and their corresponding conditional clade probability.
  void register_leafset(std::string);
  long int set2id(std::set<int> leaf_set);                                           //If the set exists, returns the set id, otherwise creates a new set id for this set and returns it.  

  //numeric
  scalar_type Bi(int n2);                                                            //Returns the total number of binary tree topologies possible given a fixed bipartition between n2 leaves on one side and Gamma_size-n2 leaves on the other side.
  scalar_type Tri(int n2,int n3);                                                    //Returns the total number of binary tree topologies possible given a fixed trifurcation between n2 leaves in one clade, n3 in another, and Gamma_size-n2-n3 leaves in the last one.

  scalar_type p_dip(long int g_id,long int gp_id,long int gpp_id);                   //Probability of a trifurcation given by the ids of the clades.
  scalar_type p_dip(std::set<int> gamma,std::set<int> gammap,std::set<int> gammapp); //Probability of a trifurcation given by the leaf sets of the clades.
  scalar_type p_bip(long int g_id);                                                  //Probability of a bipartition given by its id. Uses the correction term alpha.
  scalar_type p_bip(std::set<int> gamma);                                            //Probability of a bipartition given by its leaf set. 

  //nuisance
  std::string set2name(std::set<int> leaf_set);                                      //Prints the leaf names of leaves contained in a leaf set
  std::string random_split(std::set<int> gamma);                                     //Recursive function that returns a random subtree given a leaf set as input and given the approx_posterior object. Can return clades never observed in the posterior sample.
  std::vector<std::string>  all_trees(std::set <int > gamma);//del-loc
  scalar_type count_trees();
  scalar_type count_trees(long int g_id);
  scalar_type count_all_trees(std::set <int> gamma);
	void setAlpha ( scalar_type a );
	void setBeta ( scalar_type b );

 private:
  ;
};



//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
//"Given a set S, the power set (or powerset) of S, written P(S), or 2S, is the set of all subsets of S."
template<typename Set> std::set<Set> powerset(const Set& s, size_t n);
//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
template<typename Set> std::set<Set> powerset(const Set& s);
