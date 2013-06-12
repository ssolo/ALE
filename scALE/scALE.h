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
	std::vector<double> NeTs; 									//branch-wise parameters of the coalescent process, one per branch.
	
	std::map < long int, std::string > geneCladeIdToSpecies;                                    //Map between clade id (from the gene approx_posterior object) and species included in that clade.
	std::map < int, std::vector < map < pair <int, int  >, int > > > branch_counts;             //del-loc

	//std::map < long int, std::pair< long int, std::pair < scalar_type > > > q; //del-loc. Map between resolved clade id (from the species tree approx_posterior object), and pair between clade id from the gene tree approx_posterior object and a vector containing the probability of observing the gene tree clade at each slice of the species tree branch according to the scALE model.
	std::vector < std::pair < long int, std::vector < scalar_type > > > q; //del-loc. Same as above, but instead of a map we use a vector.	
	
	
	void computeProbabilityOfCladeInSpeciesTreeBranch (int gCladeId, 
													   long int speciesTreeResolution, 
													   int numberOfSlicesPerBranch,
													   std::vector < std::pair < long int, scalar_type > > q ) ; //Computes the probability of a given gene tree clade in a given species tree branch, at all time slices except 0.
	
	void computeProbabilityOfCladeAtBeginningOfSpeciesTreeBranch (int gCladeId, 
																		 long int speciesTreeResolution,
																		 int numberOfSlicesPerBranch, 
																  std::vector < std::pair < long int, std::vector < scalar_type > > > q ) ; //Computes the probability of a given gene tree clade at the beginning of a given species tree non-leaf branch, at time slice 0.


	
public:
	void construct(); //Constructs an object given a species tree and population size.
	scALE();
	~scALE();
	scalar_type p(approx_posterior *gale);                            //Computes the probability of an approx_posterior corresponding to a gene tree distribution according to the species tree distribution and parameter values.

	
}

#endif
