#include "ALE.h"
#include "ALE_util.h"
using namespace std;

scalar_type count_trees(approx_posterior * ale)
{
  //For clade g, a split is defined by the the complementary subclades gp and gpp
  //..similar to the ALE manuscripts (gamma',gamma''|gamma) notation.. 
  //..or put differently C[1]==g,C[2]==gp and C[3]==gp ..      

  // The loop below goes over all clades, i.e. the same as line 3 in algo 2 ( ..in increasing size order do:) 
  // g_id_count stores the number of trees (subtrees of G) that can be amalgamated for the clade g.
  // In the approx_posterior object each clade has an id, a long int, that I call, e.g. g_id 
  map<long int,scalar_type> g_id_count;//del-loc
  // ale->size_ordered_bips is a map: int "size of clade" -> vector <long int> "vector of clade ids"
  for (map <int, vector <long int > > :: iterator it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	long int g_id=(*jt);

	// leaves have size one, and only one tree can be amalgamated
	if ((*it).first==1)
	  g_id_count[g_id]=1;
	else
	  {
	    g_id_count[g_id]=0;
	    // for non-leaves we have the recursion 
	    // n(g) = sum "over splits of g" n(gp)*n(gpp) 
	    // the ale->Dip_counts object is a bit confusing, but it was the most efficient way I could find to store splits  
	    // ale->Dip_counts is a map: long int "clade id" -> map < set<long int>,scalar_type> "a map object recording a split and the number of time we saw it".  
	    // the map < set<long int>,scalar_type> "a map object recording a split and the number of time we saw it"
	    // consist of a set which has always two parts, gp and gpp (i.e. C[2] and C[3]) and the scalar_type is a double recording the number of times we saw this split 
	    // in short the loop below implements << sum "over splits of g" >>
	    for (map< set<long int>,scalar_type> :: iterator kt = ale->Dip_counts[g_id].begin(); kt != ale->Dip_counts[g_id].end(); kt++)
	      {	  
		vector <long int> parts;
		for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
		long int gp_id=parts[0];
		long int gpp_id=parts[1];	    
		// n(g) += n(gp)*n(gpp) 
		g_id_count[g_id]+=g_id_count[gp_id]*g_id_count[gpp_id];
	      }
	  }
      }

  scalar_type count_trees_g=0;  
  // I handle separately splits that correspond to roots of G
  // I record the number of times I saw a clade such that its father was the root.. 
  // the loop below implements << sum "over splits of g" >> for splits corresponding to roots of G  
  // ale->Bip_counts map: long int "clade id" -> scalar_type "number of times we saw this split" 
  // here this is like C[2] -> number of times C[2] was present such that C[1] = all genes..  
  // with the condition that C[2] is the smaller subclade of L(G), else if |C[2]|=|C[3]| we record both C[2] and C[3] taking care to consider this below (messy sorry)   
  for (map <long int,scalar_type> :: iterator it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
    {
      long int g_id=(*it).first;
      // ale->id_sets is a map: long int "clade id" -> set<int> "set of gene ids"
      set <int> gamma=ale->id_sets[g_id];
      // since this is a root L(G)/C[2]=C[3]
      set <int> not_gamma;
      for (set<int>::iterator st=ale->Gamma.begin();st!=ale->Gamma.end();st++)
	if (gamma.count(*st)==0)
	  not_gamma.insert(*st);
      
      if ( gamma.size()>not_gamma.size())
	// n(g) += n(gp)*n(gpp) 
	count_trees_g+=g_id_count[ale->set_ids[gamma]]*g_id_count[ale->set_ids[not_gamma]];
      else if (gamma.size()==not_gamma.size())
	// for the case where |C[2]|=|C[3]|, we divide by two to compensate for over-counting 
	count_trees_g+=g_id_count[ale->set_ids[gamma]]*g_id_count[ale->set_ids[not_gamma]]/2.0;
    }
  g_id_count.clear();
  return count_trees_g;
}


int main(int argc, char ** argv)
{

  string ale_file=argv[1];
  string ale_name=ale_file+".ale";
  approx_posterior * ale;
  int burnin=0;
  if (argc>2) 
    burnin=atoi(argv[2]);
  ale=observe_ALE_from_file(ale_file,burnin);
  cout << "# observe "<< ale->observations << "trees from: " <<  argv[1] << endl;
  ale->save_state(ale_name);
  cout << "# saved in "<< ale_name<<endl;

  // we delete it and load it back below
  //some info about our ale:
  cout << "Read summary of tree sample for "<<ale->observations<<" trees from: " << ale_name <<".. with :" << ale->count_trees() << " possible amalgamations .." << endl << endl ;

  delete ale;
  
  ale=load_ALE_from_file(ale_name);
  
  // above I put a version of the count_trees function which has comments to explain the interface the approx_posterior class provides for iterating over splits, it should give the same result as the previous count

  cout << " counting again ..  :" << count_trees(ale) << " possible amalgamations .." << endl << endl ;
  return 1;

}
