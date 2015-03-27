//
//  scALE.h
//  ALEODTLSpeciesTree
//
//  Created by Bastien Boussau on 07/06/13.
//  Copyright 2013 UC Berkeley. All rights reserved.
//

#ifndef scALE_h
#define scALE_h

#include "ALE.h"

/****************************************************************************
 // scALE class.
 // This class provides functions for computing the probability
 // of a distribution of species trees, given a distribution of gene trees and
 // branchwise parameters of the coalescent.
 // Tree distributions are summarized by conditional clade probabilities.
 // Using Dynamic programming, we traverse both the species tree distribution 
 // and the gene tree distribution, one pair of clades at a time.
 // This class contains an approx_posterior object for the species tree distribution,
 // with branch-wise parameters of population size Ne times time T (one parameter NeT per branch),  
 // and an approx_posterior object for the gene tree distribution. 
 *****************************************************************************/
class exODT_model
{
private:
	
	std::map <std::string,scalar_type> scalar_parameter;//del_loc
	std::map <std::string,std::vector <scalar_type> > vector_parameter;//del_loc
	std::map <std::string,std::string> string_parameter;//del_loc
	
	std::vector <std::string> speciesNames;

	
	approx_posterior * sale_pointer;                            //Pointer to an approx_posterior object used to describe a species tree distribution. Used for dynamic programming in p for instance.	
	approx_posterior * gale_pointer;                            //Pointer to an approx_posterior object used to describe a gene tree distribution. Used for dynamic programming in p for instance.
	std::map < std::set < long int >, double> NeTs; 									//branch-wise parameters of the coalescent process, one per resolved clade of the species tree.
	
	std::map < long int, std::string > geneCladeIdToSpecies;                                    //Map between clade id (from the gene approx_posterior object) and species included in that clade.
	std::map < int, std::vector < map < pair <int, int  >, int > > > branch_counts;             //del-loc

	//std::map < long int, std::pair< long int, std::pair < scalar_type > > > q; //del-loc. Map between resolution of a clade (from the species tree approx_posterior object), and pair between clade id from the gene tree approx_posterior object and a vector containing the probability of observing the gene tree clade at each slice of the species tree branch according to the scALE model.
	std::map < std::set< long int >,  std::pair < long int, std::vector < scalar_type > > > q; //del-loc. Same as above, but instead of a map we use a vector.	
	
	/******************************************************************************
	 //Computes the probability of a given gene tree clade in a given species tree branch, at all time slices except 0.
	 //Uses formula: \frac{d P(\gamma, s, t)}{dt} = -\theta_s \sum_{\bar g} \Pi_{\bar g} \sum_{s',s''} \Pi_{s'} \Pi_{s''} P(\bar \gamma, s', t) P(\bar \gamma, s'', t) \\ + \theta_s \sum_{ g} \Pi_{ g} \sum_{s',s''}\Pi_{s'}\Pi_{s''}P( \gamma', s', t)P(\gamma'', s'', t)
	 //More clearly: the probability of seeing a gene tree clade gamma at time slice n+1 is equal to the probability at time slice n of seeing two daughter clades which then coalesce into gamma before time slice n minus the probability that gamma has coalesced with another clade, and thus disappeared.
	 *******************************************************************************/
	void computeProbabilityOfCladeInSpeciesTreeBranch (int gCladeId, 
													   std::set < long int > speciesTreeResolution, 
													   int numberOfSlicesPerBranch) ; 
	
	/******************************************************************************
	//Computes the probability of a given gene tree clade at the beginning of a given species tree non-leaf branch, at time slice 0.
	//Uses formula: P(\gamma, s, t_{bottom}^s)=\sum_{s'} \Pi_{s'} P (\gamma, s', t_{top}^{s'})
	*******************************************************************************/
	void computeProbabilityOfCladeAtBeginningOfSpeciesTreeBranch (int gCladeId, 
																  std::set < long int > speciesTreeResolution,
																  int numberOfSlicesPerBranch) ; 



	
public:
	void construct(); //Constructs an object given a species tree and population size.
	scALE();
	~scALE();
	scalar_type p(approx_posterior *gale);                            //Computes the probability of an approx_posterior corresponding to a gene tree distribution according to the species tree distribution and parameter values.

	
}

#endif
