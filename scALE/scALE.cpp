//
//  scALE.cpp
//  Coale
//
//  Created by Bastien Boussau on 07/06/13.
//  Copyright 2013 UC Berkeley. All rights reserved.
//

#include <iostream>
#include "scALE.h"

scALE::scALE()
{
	//some default parameters
	string_parameter["gene_name_separators"]="_@";
	scalar_parameter["species_field"]=0;
	scalar_parameter["event_node"]=0;
	scalar_parameter["min_bip_count"]=-1;
	scalar_parameter["min_branch_lenghts"]=0;  
	// length of "stem" branch above root
	scalar_parameter["stem_length"]=1;
	//number of subdiscretizations for ODE calculations
	//Corresponds to maximum number of coalescences on a given branch of the species tree.
	scalar_parameter["DD"]=10;
}


void scALE::construct( approx_posterior *sale )
{
	sale_pointer = sale;
	speciesNames = sale_pointer->getLeafNames();
	
}





scalar_type scALE::p(approx_posterior *gale) {
	gale_pointer = gale;
	approx_posterior *sale = sale_pointer;
	//directed partitions and their sizes, for the gene and the species tree distributions 
	vector <long int>  g_ids;       //del-loc. Vector of leaf set (=clade) ids for the 
									//gene tree distribution, ordered by their size, small to large.
	vector <long int>  g_id_sizes;  //del-loc. Numbers of leaves in the above sets. 
	
	vector <long int>  s_ids;       //del-loc. Vector of leaf set (=clade) ids for the 
									//species tree distribution, ordered by their size, small to large.
	vector <long int>  s_id_sizes;  //del-loc. Numbers of leaves in the above sets. 

	gale->computeOrderedVectorOfClades ( g_ids, g_id_sizes );	
	sale->computeOrderedVectorOfClades ( s_ids, s_id_sizes );	
	size_t numSpeciesClades = s_ids.size();
	size_t numGeneClades = g_ids.size();
	size_t numberOfSlicesPerBranch = scalar_parameter["DD"];

	
	//Need to empty q
	for ( std::vector < std::pair < long int, std::vector < scalar_type > > > q::iterator it=q.begin() ; it!=q.end() ; ++it )
    {
		for (size_t slice=0; slice < numberOfSlicesPerBranch ; ++slice) {
			(*it).second[i] = 0.0;
		}
    }      
	if ( q.size () != numSpeciesClades )
		q.resize ( numSpeciesClades );
	
	
	// gene<->species mapping
	for ( int i=0; i<(int)g_ids.size(); i++ ) //Going through each clade of the gene approx_posterior
    {
		long int g_id=g_ids[i];				
		if (g_id_sizes[i]==1) //a leaf, mapping is by name
		{
			string gene_name = gale->id_leaves[(* (gale->id_sets[g_id].begin()) )];
			vector <string> tokens;
			boost::split(tokens,gene_name,boost::is_any_of(string_parameter["gene_name_separators"]),boost::token_compress_on);
			string species_name;
			if ((int)scalar_parameter["species_field"]==-1)
				species_name=tokens[tokens.size()-1];	  
			else
				species_name=tokens[(int)scalar_parameter["species_field"]];	  
			geneCladeIdToSpecies[g_id]=species_name;
		}	
		else {
			break;
		}
    }
	

	//Now, the main loop, iterating over the clades of the species tree distribution
	long int spCladeId, gCladeId;
	map< set<long int>,scalar_type> speciesCladeResolutions ;
	map< set<long int>,scalar_type> geneCladeResolutions ;
	for ( size_t i = 0 ; i < numSpeciesClades ; ++i ) 
	{
		spCladeId = s_ids[i];
		speciesCladeResolutions = sale->Dip_counts[spCladeId];
		//Second loop, over the resolutions of the species tree clade spCladeId.
		for (map< set<long int>,scalar_type> :: iterator spResolution = speciesCladeResolutions.begin(); spResolution != speciesCladeResolutions.end(); ++spResolution) //Going through all resolutions of the clade spCladeId
		{	 
			//Third loop, over the clades of the gene tree distribution.
			for ( size_t j = 0 ; j < numGeneClades ; ++j ) 
			{
				gCladeId = 	g_ids[j];
				geneCladeResolutions = gale->Dip_counts[gCladeId];
				//Fourth loop, over the resolutions of the gene tree clade gCladeId.
				for (map< set<long int>,scalar_type> :: iterator gResolution = geneCladeResolutions.begin(); gResolution != geneCladeResolutions.end(); ++gResolution) //Going through all resolutions of the clade gCladeId
				{
					if (s_id_sizes[i] == 1) 
					{//The species tree clade is a leaf
						std::string currentSpeciesId = sale->id_leaves[(* (sale->id_sets[spCladeId].begin()) )];
						if (g_id_sizes[i] == 1) 
						{//The gene tree clade is a leaf
							if (sale->id_leaves[(* (sale->id_sets[spCladeId].begin()) )] == geneCladeIdToSpecies[gCladeId] )
							{//Gene corresponds to the species
								q[ *(spResolution)->first ][gCladeId][0] = 1.0;
								computeProbabilityOfCladeInSpeciesTreeBranch (gCladeId, *(spResolution)->first,numberOfSlicesPerBranch, q );
							}
							else 
							{
								for (size_t slice=0; slice < numberOfSlicesPerBranch ; ++slice) {
									q[ *(spResolution)->first ][gCladeId][slice] = 0.0;	
								}
							}
						}//End the gene tree clade is a leaf
						else 
						{//The gene tree clade is not a leaf
							//Check that all leaves in the gene tree clade correspond to the species
							std::set <int> geneLeafIds = gale->id_sets[gCladeId];
							bool allGenesFromCurrentSpecies = true;
							for (std::set<int>::iterator it =  geneLeafIds.begin(); it!= geneLeafIds.end() ++it ) {
								if (geneCladeIdToSpecies [*it] != currentSpeciesId ) {
									allGenesFromCurrentSpecies=false;
									break;
								}
							}
							if (allGenesFromCurrentSpecies) {
								//need to compute the probability of observing the given clade in this species tree branch
								q[ *(spResolution)->first ][gCladeId][0] = 0.0;	//We can't have an entire clade already coalesced at sampling time!
								computeProbabilityOfCladeInSpeciesTreeBranch (gCladeId, *(spResolution)->first, numberOfSlicesPerBranch, q );
							}
							else {
								for (size_t slice=0; slice < numberOfSlicesPerBranch ; ++slice) {
									q[ *(spResolution)->first ][gCladeId][slice] = 0.0;
								}
							}
						}//End the gene tree clade is not a leaf
					}//End the species tree clade is a leaf
					else 
					{ //The species tree clade is not a leaf
						std::set< std::string > speciesInClade;
						std::set< int > speciesLeafIds = sale->id_sets[spCladeId];
						for (std::set<int>::iterator it =  speciesLeafIds.begin(); it!= speciesLeafIds.end() ++it ) {
							speciesInClade.insert ( sale->id_leaves[(*it)] );
						}
						if (g_id_sizes[i] == 1) 
						{//The gene tree clade is a leaf
							computeProbabilityOfCladeAtBeginningOfSpeciesTreeBranch (gCladeId, *(spResolution)->first, numberOfSlicesPerBranch, q );
							if (q[ *(spResolution)->first ][gCladeId][0] == 0) { //The current gene  comes from a species not in this clade
								for (size_t slice = 1; slice < numberOfSlicesPerBranch ; ++slice) {
									q[ *(spResolution)->first ][gCladeId][slice] = 0.0;
								}
							}
							else {
								computeProbabilityOfCladeInSpeciesTreeBranch (gCladeId, *(spResolution)->first, numberOfSlicesPerBranch, q );
							}
						}//End the gene tree clade is a leaf
						else 
						{//The gene tree clade is not a leaf
							//Check that all leaves in the gene tree clade correspond to the species in the species tree clade
							std::set <int> geneLeafIds = gale->id_sets[gCladeId];
							bool allGenesFromCurrentSpecies = true;
							for (std::set<int>::iterator it =  geneLeafIds.begin(); it!= geneLeafIds.end() ++it ) {
								if (speciesInClade.find ( geneCladeIdToSpecies [*it] ) == speciesInClade.end() ) {
									allGenesFromCurrentSpecies=false;
									break;
								}
							}
							if (allGenesFromCurrentSpecies) {
								computeProbabilityOfCladeAtBeginningOfSpeciesTreeBranch (gCladeId, *(spResolution)->first, numberOfSlicesPerBranch, q );
								computeProbabilityOfCladeInSpeciesTreeBranch (gCladeId, *(spResolution)->first, numberOfSlicesPerBranch, q );
							}
							else {
								for (size_t slice=0; slice < numberOfSlicesPerBranch ; ++slice) {
									q[ *(spResolution)->first ][gCladeId][slice] = 0.0;
								}
							}
						}//End the gene tree clade is not a leaf
					}//End the species tree clade is not a leaf
				}//End loop over gene tree clade resolutions
			}//End loop over gene tree clades
		}//End loop over species tree clade resolutions
	}//End loop over species tree clades
	
	//del-locs
	g_ids.clear();
	g_id_sizes.clear();
	
	return root_sum;	
}



void scALE::computeProbabilityOfCladeInSpeciesTreeBranch (int gCladeId, 
															std::set < long int >  speciesTreeResolution,
															int numberOfSlicesPerBranch, 
															std::vector < std::pair < long int, std::vector < scalar_type > > > q ) {
	for (size_t slice=1; slice < numberOfSlicesPerBranch ; ++slice) {
		//TODO !
	}
	return;
}



void scALE::computeProbabilityOfCladeAtBeginningOfSpeciesTreeBranch (int gCladeId, 
																	 std::set < long int > speciesTreeResolution,
																	 int numberOfSlicesPerBranch, 
																	 std::vector < std::pair < long int, std::vector < scalar_type > > > q ) {
	q[ speciesTreeResolution ][gCladeId][0] = 0.0;
	int lastSlice = numberOfSlicesPerBranch - 1;
	for (std::set<long int>::iterator speciesTreeDaughterClade = speciesTreeResolution.begin() ; speciesTreeDaughterClade != speciesTreeResolution.end() ; ++speciesTreeDaughterClade) {
		speciesTreeDaughterCladeResolutions = sale->Dip_counts[*(speciesTreeDaughterClade)];
		//Second loop, over the resolutions of the species tree clade speciesTreeDaughterClade.
		for (map< set<long int>,scalar_type> :: iterator spDaughterResolution = speciesTreeDaughterCladeResolutions.begin(); spDaughterResolution != speciesTreeDaughterCladeResolutions.end(); ++spDaughterResolution) //Going through all resolutions of the clade speciesTreeDaughterClade
		{
			q[ speciesTreeResolution ][gCladeId][0] += *(spDaughterResolution)->second * q[ *(spDaughterResolution)->first ][gCladeId][lastSlice];
		}
	}
	return;
}
