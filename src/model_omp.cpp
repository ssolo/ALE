#include "exODT.h"
#include "omp.h"
using namespace std;
using namespace bpp;

//#include <bitset>

//oMP// 
//oMP// openMP added into this function
//oMP// 

scalar_type exODT_model::p(approx_posterior *ale)
{
  ale_pointer=ale;
  //directed partitions and their sizes
  vector <long int>  g_ids;       //del-loc. Vector of leaf set (=clade) ids, ordered by their size, small to large.
  vector <long int>  g_id_sizes;  //del-loc. Numbers of leaves in the above sets. 

   //First, cleaning q. 
  for (std::map<long int, std::map< scalar_type, std::map<int, scalar_type> > >::iterator it=q.begin();it!=q.end();it++)
    {
      for ( std::map< scalar_type, std::map<int, scalar_type> >::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	(*jt).second.clear();
      (*it).second.clear();
    }      
  q.clear();

  //cout << "start" << endl;  
  //iterate over directed partitions (i.e. clades) ordered by the number of leaves
  //cout << "start loop" << endl;
  //test  
  //long int tmp_g_id=-1;
  //cout << ale->set2name(ale->id_sets[tmp_g_id]) <<endl;
  //test  

  //oMP// 
  //oMP// I sort the directed partitions by size (number of gene tree leaves) to ensure that we calculate things in the proper order (smaller to larger)
  //oMP// 

  for (map <int, vector <long int > > :: iterator it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	g_ids.push_back((*jt));
	g_id_sizes.push_back((*it).first);
      }
  //root bipartition needs to be handled separately (and last, given it's the largest)
  g_ids.push_back(-1);
  g_id_sizes.push_back(ale->Gamma_size);

  /*
  // gene<->species mapping
  for (int i=0;i<(int)g_ids.size();i++) //Going through each clade of the approx_posterior
    {
      long int g_id=g_ids[i];
      cerr<<"g_id: "<<g_id<<"\n";
      for (int rank=0;rank<last_rank;rank++) //Going through time slices, from leaves to root
        {
	  cerr<<"\trank: "<<rank<<"\n";
          int n=time_slices[rank].size(); //Number of branches in that time slice
          for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++) //Going through the subslices
            {
	      cerr<<"\t\tt_i: "<<t_i<<"\n";
              scalar_type t=time_slice_times[rank][t_i];
              for (int branch_i=0;branch_i<n;branch_i++)
                {
                  int e = time_slices[rank][branch_i];
		  cerr<<"\t\t\te: "<<e<<" branch_i: "<<branch_i<<"\n";
                  q[g_id][t][e]=0; //initializing q with 0s for each clade of the ale, each time subslice, each branch
                }
              q[g_id][t][alpha]=0; //I don't understand the purpose of that: alpha=-1, so it's already 0, right?
            }
        }
      
      if (g_id_sizes[i]==1) //a leaf, mapping is by name
        {
          string gene_name=ale->id_leaves[(* (ale->id_sets[g_id].begin()) )];
          vector <string> tokens;
          boost::split(tokens,gene_name,boost::is_any_of(string_parameter["gene_name_separators"]),boost::token_compress_on);
          string species_name;
          if ((int)scalar_parameter["species_field"]==-1)
            species_name=tokens[tokens.size()-1];
          else
            species_name=tokens[(int)scalar_parameter["species_field"]];
          gid_sps[g_id]=species_name;
        }
    }
  */

  vector<vector<vector<map<int, scalar_type> > > > qvec;
  
  // gene<->species mapping
  //  for (int i=0;i<(int)g_ids.size();i++) //Going through each clade of the approx_posterior
  // {
  //    long int g_id=g_ids[i];
  //   cerr<<"i: "<<i<<" g_id: "<<g_id<<"\n";
  for (int g_id= -1; g_id<(int)g_ids.size(); g_id++) // ACCES par g_id+1
    {
      //cerr<<"g_id: "<<g_id<<"\n";
      if (g_id==0){ // n'existe pas --> case vide
	vector<vector<map<int, scalar_type> > > vrank;
	vector<map<int, scalar_type> > vt_i;
	map<int, scalar_type> vbranch;
	vt_i.push_back(vbranch);
	vrank.push_back(vt_i);
	qvec.push_back(vrank);
      }
      else{
	//vector<vector<vector<scalar_type> > > vrank;
	vector<vector<map<int, scalar_type> > > vrank;
	for (int rank=0;rank<last_rank;rank++) //Going through time slices, from leaves to root
	  {
	    //cerr<<"\trank: "<<rank<<"\n";
	    int n=time_slices[rank].size(); //Number of branches in that time slice
	    //vector<vector<scalar_type> > vt_i;
	    vector<map<int, scalar_type> > vt_i;
	    for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++) //Going through the subslices
	      {
		//cerr<<"\t\tt_i: "<<t_i<<"\n";
		//scalar_type t=time_slice_times[rank][t_i];
		//vector<scalar_type> vbranch(n, 0.);
		map<int, scalar_type> vbranch;
		for (int branch_i=0;branch_i<n;branch_i++)
		  {	    
		    int e = time_slices[rank][branch_i];
		    //cerr<<"\t\t\te: "<<e<<" branch_i: "<<branch_i<<"\n";
		    
		    //q[g_id][t][e]=0; //initializing q with 0s for each clade of the ale, each time subslice, each branch
		    vbranch[e]=0;
		  }
		//q[g_id][t][alpha]=0; //I don't understand the purpose of that: alpha=-1, so it's already 0, right? 
		vbranch[alpha]=0;
		vt_i.push_back(vbranch);
	      }
	    vrank.push_back(vt_i);
	  }
	qvec.push_back(vrank); // g_id+1 to manage -1
      }
    }  
  for (int i=0;i<(int)g_ids.size();i++) //Going through each clade of the approx_posterior
    {
      long int g_id=g_ids[i];
      if (g_id_sizes[i]==1) //a leaf, mapping is by name
	{
      /*  int id = 0;
         boost::dynamic_bitset<> temp = ale->id_sets[g_id];
        for (auto i = 0; i < ale->Gamma_size + 1 ; ++i) {
//            if ( BipartitionTools::testBit ( temp, i) ) {
            if ( temp[i] ) {
                id = i;
                break;
            }
        }*/
        
        int id = 0;
        for (auto i=0; i< ale->Gamma_size + 1; ++i) {
            if ( ale->id_sets[g_id][i] ) {
                id=i;
                break;
            }
        }
        
        string gene_name=ale->id_leaves[ id /*g_id*/ ];

//        string gene_name=ale->id_leaves[ g_id ];
	//  string gene_name=ale->id_leaves[(* (ale->id_sets[g_id].begin()) )];
	  vector <string> tokens;
	  boost::split(tokens,gene_name,boost::is_any_of(string_parameter["gene_name_separators"]),boost::token_compress_on);
	  string species_name;
	  if ((int)scalar_parameter["species_field"]==-1)
	    species_name=tokens[tokens.size()-1];	  
	  else
	    species_name=tokens[(int)scalar_parameter["species_field"]];	  
	  gid_sps[g_id]=species_name;
	}	 
    }  
   
  //oMP// 
  //oMP// below is the loop that iterates over the sorted g_ids, it is this one that should be amicable to openMP  
  //oMP// the important thing is that we can only do the g_ids in parallel that have the same number of leaves
  //oMP// hence the sorting above..
  //oMP// 
  //oMP// the calculation fills out the global q, cf. exODT.h, this is later needed for sampling reconciliations! 
  //oMP// 
  //oMP// 


  std::map<int, vector<int> > size2i; //Map between clade size and vector of ids of the clades of that size. Seems to me this size2i could be built once and for all, and saved in the approx_posterior object, and not reconstructed every time?
    
  for (int i=0;i<(int)g_ids.size();i++) { //Filling up size2i
    if (size2i.count(g_id_sizes[i]) == 0)
      size2i[g_id_sizes[i]] = vector<int> ();                  
    size2i[g_id_sizes[i]].push_back(i);
  }
        
  for (map <int, vector < int > > :: iterator it2 = size2i.begin(); it2 != size2i.end(); it2++) 
    { //
      int j=0;
      int siz = (int)it2->second.size(); //Number of clades with that size
      if (siz <= 4 )// num_threads ) //If few clades: inside loop parallelization
        {
	  for ( j=0 ; j < siz ;j++)
	    {
	      int i = it2->second.at(j);
             
	      // directed partition (dip) gamma's id  
	      bool is_a_leaf=false;
	      long int g_id=g_ids[i];	
	      if (g_id_sizes[i]==1)
		is_a_leaf=true;
             
	      vector <long int> gp_ids;//del-loc
	      vector <long int> gpp_ids;//del-loc
	      vector <scalar_type> p_part;//del-loc
	      if (g_id!=-1)
		{
		  {
		    for (unordered_map< pair<long int, long int>,scalar_type> :: iterator kt = ale->Dip_counts[g_id].begin(); kt != ale->Dip_counts[g_id].end(); kt++)
		      {	  
			/*			vector <long int> parts;
			for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
			long int gp_id=parts[0];
			long int gpp_id=parts[1];	    */
			pair<long int, long int> parts = (*kt).first; 
			long int gp_id = parts.first;
			long int gpp_id = parts.second;
			gp_ids.push_back(gp_id);
			gpp_ids.push_back(gpp_id);
			if (ale->Bip_counts[g_id]<=scalar_parameter["min_bip_count"])
			  p_part.push_back(0);
			else
			  p_part.push_back(ale->p_dip(g_id,gp_id,gpp_id));
		      }
		  }
		}
	      else
		{
		  //root bipartition needs to be handled separately
		  map<set<long int>,int> bip_parts;
		  for (map <long int,scalar_type> :: iterator it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
		    {
		      long int gp_id=(*it).first;
		      boost::dynamic_bitset<> gamma = ale->id_sets.at(gp_id);
		      boost::dynamic_bitset<> not_gamma = ~gamma;
                not_gamma[0] = 0;
		      /*for (set<int>::iterator st=ale->Gamma.begin();st!=ale->Gamma.end();st++)
			if (gamma.count(*st)==0)
			  not_gamma.insert(*st);*/
		      long int gpp_id = ale->set_ids.at(not_gamma);
		      set <long int> parts;
		      parts.insert(gp_id);
		      parts.insert(gpp_id);
		      bip_parts[parts]=1;
		      /*gamma.clear();
		      not_gamma.clear();*/
		    }
		  for (map<set<long int>,int> :: iterator kt = bip_parts.begin();kt!=bip_parts.end();kt++)
		    {
		      vector <long int> parts;
		      for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
		      long int gp_id=parts[0];
		      long int gpp_id=parts[1];	    
		      gp_ids.push_back(gp_id);
		      gpp_ids.push_back(gpp_id);
                //Here we can create a new ale->Bip_counts[gp_id], in particular for leaves.
                //We may want to add the leaf entries for Bip_counts when Bip_counts is first created.
		      if (ale->Bip_counts[gp_id]<=scalar_parameter.at("min_bip_count") and not ale->Gamma_size<4)
			p_part.push_back(0);	      
		      else
			p_part.push_back(ale->p_bip(gp_id));	      
		    }
		  bip_parts.clear();
		}
             
	      int N_parts=gp_ids.size();
                 
	      //iterate over all positions along S
	      for (int rank=0;rank<last_rank;rank++)
		{
		  int n=time_slices[rank].size();
		  for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
		    {
                         
		      //######################################################################################################################
		      //#########################################INNER LOOP##################################################################
		      //######################################################################################################################
                     
		      scalar_type t=time_slice_times[rank][t_i];
		      scalar_type tpdt;//,tpdt_nl;
		      int tpdt_rank, tpdt_t_i;
		      if ( t_i <  (int)time_slice_times[rank].size()-1 ){
			tpdt=time_slice_times[rank][t_i+1];
			tpdt_rank = rank;
			tpdt_t_i = t_i+1;
		      }
		      else if (rank<last_rank-1){
			tpdt=time_slice_times[rank+1][0];
			tpdt_rank = rank+1;
			tpdt_t_i = 0;
		      }
		      else{
			//top of root stem
			tpdt=t_begin[time_slices[rank][0]]; //PBM PBM PBM
			tpdt_rank = rank;  //PBM PBM PBM
			tpdt_t_i = 0;	 //PBM PBM PBM		
		      }
                     
		      bool tpdt_nl_is_t = true;
		      if (scalar_parameter["event_node"]==1 and false)
			;//tpdt_nl=t;			
		      else{
			//tpdt_nl=tpdt;
			tpdt_nl_is_t = false;
		      }
                     
		      //root
		      scalar_type Delta_t=tpdt-t;
		      //scalar_type N=vector_parameter["N"][rank];
		      scalar_type Delta_bar=vector_parameter["Delta_bar"][rank];
		      //scalar_type Lambda_bar=vector_parameter["Lambda_bar"][rank];
		      //OMG
		      //scalar_type p_Delta_bar=1-exp(-Delta_bar/N*Delta_t);			     
		      scalar_type p_Delta_bar=Delta_bar*Delta_t;			     
		      scalar_type Ebar=Ee[-1][t];
                     
		      //boundaries for branch alpha virtual branch  
		      //boundary at present
		      if (t==0) 
			{
			  //#pragma omp critical 
			  {
			    //q[g_id][t][alpha]=0;
			    qvec[g_id+1][rank][t_i][alpha]=0;
			  }
			}
		      //boundary between slice rank and rank-1 slice is trivial	
		      ;//q[g_id][t][alpha]=q[g_id][t][alpha];	  
		      //boundaries for branch alpha virtual branch.  
		      if(1)
			{
#pragma omp parallel for schedule(dynamic,1) //p2
			  for (int branch_i=0;branch_i<n;branch_i++)
			    {	  
			      {
				int e = time_slices[rank][branch_i];
                               
				//boundaries for branch e
				//boundary at present
				if (t==0)
				  {
				    if (is_a_leaf && extant_species[e]==gid_sps[g_id]) {	
				      //#pragma omp critical 
				      {
					//q[g_id][t][e]=1;
					qvec[g_id+1][rank][t_i][e]=1;
				      }
				    }
				    else {
				      //#pragma omp critical 
				      {			
					//[i][rank][t_i][branch_i]=0;		
					//q[g_id][t][e]=0;
					qvec[g_id+1][rank][t_i][e]=0;
				      }
				    }
				  }
				//boundary between slice rank and rank-1
				else if (t_i==0)
				  {
				    //terminating branch is last in time_slices and defines a represented speciation 
				    if (branch_i==n-1 && rank>0)
				      {
					int f=daughters[e][0];
					int g=daughters[e][1];

					scalar_type Eft=Ee[f][t];
					scalar_type Egt=Ee[g][t];
                                         
					scalar_type q_sum=0;
					//q[g_id][t][e]=0;
                                         
					scalar_type SL_fLg=qvec[g_id+1][rank][t_i][f]*Egt;
					scalar_type SL_Lfg=qvec[g_id+1][rank][t_i][g]*Eft;
					//scalar_type SL_fLg=q[g_id][t][f]*Egt;
					//scalar_type SL_Lfg=q[g_id][t][g]*Eft;
					//SL EVENT
					q_sum+=SL_fLg+SL_Lfg;
					//q[g_id][t][e]=q[g_id][t][f]*Egt + q[g_id][t][g]*Eft;
					//SL.
                                         
					//non-leaf directed partition
					if (not is_a_leaf)
					  for (int i=0;i<N_parts;i++)
					    {	
					      long int gp_id=gp_ids[i];
					      long int gpp_id=gpp_ids[i];	    
					      scalar_type pp=p_part[i];
					      scalar_type S_pf_ppg=qvec[gp_id+1][rank][t_i][f]*qvec[gpp_id+1][rank][t_i][g]*pp;
					      scalar_type S_ppf_pg=qvec[gpp_id+1][rank][t_i][f]*qvec[gp_id+1][rank][t_i][g]*pp;
					      //scalar_type S_pf_ppg=q[gp_id][t][f]*q[gpp_id][t][g]*pp;
					      //scalar_type S_ppf_pg=q[gpp_id][t][f]*q[gp_id][t][g]*pp;
					      //S EVENT
					      //q[g_id][t][e]+=q[gp_id][t][f]*q[gpp_id][t][g] +q[gpp_id][t][f]*q[gp_id][t][g];
					      q_sum+= S_pf_ppg + S_ppf_pg;
					      //S.
                                                 
					    }
					//#pragma omp critical 
					{
					  qvec[g_id+1][rank][t_i][e]=q_sum;
					  //q[g_id][t][e]=q_sum; 
					}
                                         
				      }
				    //branches that cross to next time slice  
				    else
				      {
					//trivial
					;//q[g_id][t][e]=q[g_id][t][e];
				      }			  
				  }		   
				//boundaries for branch e.
			      }
			    }
			}
                     
		      if(1)
			{
			  //events within slice rank at time t on alpha virtual branch
			  scalar_type G_bar=Ge[-1][t];//exp(-(Delta_bar*(n-N)/N+Lambda_bar)*Delta_t );	
			  //#pragma omp critical 
			  {        
			    qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]=0;                    
			    //q[g_id][tpdt][alpha]=0;
			  }
			  scalar_type q_sum=0;
			  scalar_type q_sum_nl=0;
			  //#pragma omp parallel for schedule(dynamic,1)  reduction(+:q_sum_nl) //p3 rerunning p(ale) after sample() crashes if this active !?
			  for (int branch_i=0;branch_i<n;branch_i++)			  
			     {
			      {

				int e = time_slices[rank][branch_i];		
				scalar_type tau_e=vector_parameter["tau"][e];
				//G_bar*=exp(- tau_e*Delta_t);
                                 
				//scalar_type p_Ntau_e=1-exp(-N*tau_e*Delta_t);
				//OMG
				scalar_type p_Ntau_e=tau_e*Delta_t;
				//non-leaf directed partition
				if (not is_a_leaf)
				  for (int i=0;i<N_parts;i++)
				    {	
				      long int gp_id=gp_ids.at(i);
				      long int gpp_id=gpp_ids.at(i);	    
				      scalar_type pp=p_part.at(i);
				      scalar_type T_ep_app=p_Ntau_e*qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]*pp;
				      scalar_type T_ap_epp=p_Ntau_e*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]*pp;				    
				      //scalar_type T_ep_app=p_Ntau_e*q[gp_id][t][e]*q[gpp_id][t][alpha]*pp;
				      //scalar_type T_ap_epp=p_Ntau_e*q[gp_id][t][alpha]*q[gpp_id][t][e]*pp;
				      //T EVENT
				      q_sum_nl+=T_ep_app+T_ap_epp;
				      //q[g_id][tpdt][alpha]+=p_Ntau_e*(q[gp_id][t][e]*q[gpp_id][t][alpha]+q[gp_id][t][alpha]*q[gpp_id][t][e]);
				      //T.
				    }
			      }
			     }
			  //non-leaf directed partition
			  if (not is_a_leaf)
			    for (int i=0;i<N_parts;i++)
			      {	
				long int gp_id=gp_ids[i];
				long int gpp_id=gpp_ids[i];	    
				scalar_type pp=p_part[i];
				scalar_type Sb=p_Delta_bar*(2*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha])*pp;
				//scalar_type Sb=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha])*pp;
				//S_bar EVENT
				q_sum_nl+=Sb;
				//q[g_id][tpdt][alpha]+=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha]);
				//S_bar.
                                 
			      }	           
			  if (tpdt_nl_is_t)			    
			    qvec[g_id+1][rank][t_i][alpha]+=q_sum_nl;
			  else			    
			    qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum_nl;			    			  
			  //q[g_id][tpdt_nl][alpha]+=q_sum_nl;
#pragma omp parallel for schedule(dynamic,1)  reduction(+:q_sum) //p4
			  //#pragma omp task untied 
			  for (int branch_i=0;branch_i<n;branch_i++)			  
			    {
			      int e = time_slices[rank][branch_i];		
			      scalar_type tau_e=vector_parameter["tau"][e];
			      //OMG
			      //scalar_type p_Ntau_e=1-exp(-N*tau_e*Delta_t);
			      scalar_type p_Ntau_e=tau_e*Delta_t;
			      scalar_type TLb=p_Ntau_e*Ebar*qvec[g_id+1][rank][t_i][e];
			      //scalar_type TLb=p_Ntau_e*Ebar*q[g_id][t][e];
			      //TL_bar EVENT
			      q_sum+=TLb;
			      //q[g_id][tpdt][alpha]+=p_Ntau_e*Ebar*q[g_id][t][e];
			      //TL_bar.
			    }
			  //0 EVENT
			  scalar_type empty=G_bar*qvec[g_id+1][rank][t_i][alpha];
			  //scalar_type empty=G_bar*q[g_id][t][alpha]; 
			  q_sum+=empty;
                         
			  //q[g_id][tpdt][alpha]+=G_bar*q[g_id][t][alpha];
			  //0.
			  //max
			  /*
			    if (max_term<empty) 
			    {
			    max_term=empty;
			    }
                          */
			  //max	
			  qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum;	    
			  //q[g_id][tpdt][alpha]+=q_sum;
			  //events within slice rank at time t on alpha virtual branch.
			}
		      if(1)
			{
			  //DOES PREVENT SEGFAULT
			  for (int branch_i=0;branch_i<n;branch_i++)
			    {
			      int e = time_slices[rank][branch_i];
			      qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=0;
			      //q[g_id][tpdt][e]=0;
			    }
		       
#pragma omp parallel for schedule(dynamic,1) //p5
			  for (int branch_i=0;branch_i<n;branch_i++)
			    {	
			      int e = time_slices[rank][branch_i];
			      scalar_type Get=Ge[e][t];
			      scalar_type Eet=Ee[e][t];	
			      scalar_type delta_e=vector_parameter["delta"][e];
			      scalar_type p_delta_e=delta_e*Delta_t;
                                 
			      //events within slice rank at time t on branch e 
			      //n#pragma omp critical 
			      //n{
                                     
			      //n q[g_id][tpdt][e]=0;
			      //n}
			      scalar_type q_sum=0;
			      scalar_type q_sum_nl=0;
                                 
			      //non-leaf directed partition		   
                                 
			      if (not is_a_leaf)
				for (int i=0;i<N_parts;i++)
				  {	
				    long int gp_id=gp_ids[i];
				    long int gpp_id=gpp_ids[i];	    
				    scalar_type pp=p_part[i];
				    scalar_type qpe=qvec[gp_id+1][rank][t_i][e];
				    scalar_type qppe=qvec[gpp_id+1][rank][t_i][e];
				    //scalar_type qpe=q[gp_id][t][e];
				    //scalar_type qppe=q[gpp_id][t][e];
				    scalar_type Sb_pa_ppe= p_Delta_bar*qvec[gp_id+1][rank][t_i][alpha]*qppe*pp;
				    scalar_type Sb_pe_ppa= p_Delta_bar*qpe*qvec[gpp_id+1][rank][t_i][alpha]*pp;
				    //scalar_type Sb_pa_ppe= p_Delta_bar*q[gp_id][t][alpha]*qppe*pp;
				    //scalar_type Sb_pe_ppa= p_Delta_bar*qpe*q[gpp_id][t][alpha]*pp;
				    //S_bar EVENT
				    q_sum_nl+= Sb_pa_ppe + Sb_pe_ppa;
                                         
				    //q[g_id][tpdt][e]+=p_Delta_bar*(q[gp_id][t][alpha]*q[gpp_id][t][e]+q[gp_id][t][e]*q[gpp_id][t][alpha]);			  
				    //S_bar.
                                         
				    scalar_type D=2*p_delta_e*qpe*qppe*pp;
				    //D EVENT
				    q_sum_nl+= D;
                                         
				    //q[g_id][tpdt][e]+=p_delta_e*q[gp_id][t][e]*q[gpp_id][t][e];
				    //D.
                                         
				  }
                                 
			      scalar_type SLb=p_Delta_bar*Eet*qvec[g_id+1][rank][t_i][alpha];
			      //scalar_type SLb=p_Delta_bar*Eet*q[g_id][t][alpha];
			      //SL_bar EVENT
			      q_sum_nl+=SLb;
                                 
			      //q[g_id][tpdt][e]+=p_Delta_bar*Eet*q[g_id][t][alpha];
			      //SL_bar.
			      //#pragma omp critical 
			      {			      
				if (tpdt_nl_is_t)
				  qvec[g_id+1][rank][t_i][e]+=q_sum_nl;
				else
				  qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum_nl;				
				//q[g_id][tpdt_nl][e]+=q_sum_nl;
			      }
			      scalar_type empty=Get*qvec[g_id+1][rank][t_i][e];
			      //scalar_type empty=Get*q[g_id][t][e];
			      //0 EVENT
			      q_sum+=empty;
                                 
			      //q[g_id][tpdt][e]=Get*q[g_id][t][e];
			      //0.
			      //#pragma omp critical 
			      {
				qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum;
				//q[g_id][tpdt][e]+=q_sum;
			      }
			      //events within slice rank at time t on branch e.
                             
			    }
			}
		      //######################################################################################################################
		      //#########################################INNNER LOOP##################################################################
		      //######################################################################################################################		    
		    }
		}
             
	      //}
	      gp_ids.clear();
	      gpp_ids.clear();
	      p_part.clear();
	      //scalar_type tnow=omp_get_wtime();//t->elapsed();
	      //n cout <<  "Inside loop parallelization: Nb Partitions: " << N_parts << " iteration duration: " << (tnow-tatom) << endl; ;
	      //tatom=tnow;
            }
	}
      else { ////If more clades of a given size than number of threads: outside loop parallelization
#pragma omp parallel //num_threads(8) //p6
	{
#pragma omp for schedule(dynamic,1) //p7
	  for ( j=0 ; j<siz ;j++)
	    {
	      //   std::cout << "Num Thread: "<< omp_get_thread_num()<<std::endl;
             
	      //scalar_type tatom=omp_get_wtime(); 
             
	      //            std::cout << "HERE 3"<<std::endl;
	      //    std::cout << "j: "<<j<<std::endl;
	      //  std::cout << " and it2->first: "<<it2->first << std::endl;
             
	      // std::cout << " and : "<< it2->second.at(j) <<std::endl;
	      int i = it2->second.at(j); //working on clade i
             
	      // directed partition (dip) gamma's id  
	      bool is_a_leaf=false;
	      long int g_id=g_ids[i];	//clade i has id g_id
	      if (g_id_sizes[i]==1)
		is_a_leaf=true;
             
	      vector <long int> gp_ids;//del-loc. All clades that are left daughters of clade g_id.
	      vector <long int> gpp_ids;//del-loc. All clades that are rigt daughters of clade g_id.
	      vector <scalar_type> p_part;//del-loc. Stores all probabilities of the observed resolutions of clade g_id.
	      if (g_id!=-1) //Not at the root
		{
#pragma omp critical 
		  {
		    for (unordered_map< pair<long int, long int>,scalar_type> :: iterator kt = ale->Dip_counts[g_id].begin(); kt != ale->Dip_counts[g_id].end(); kt++) //Going through all resolutions of the clade g_id
		      {	  
			/*vector <long int> parts;
			for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
			long int gp_id=parts[0];
			long int gpp_id=parts[1];	    */
			pair<long int, long int> parts = (*kt).first;
			long int gp_id = parts.first;
			long int gpp_id = parts.second;
			gp_ids.push_back(gp_id);
			gpp_ids.push_back(gpp_id);
			if (ale->Bip_counts[g_id]<=scalar_parameter["min_bip_count"])
			  p_part.push_back(0);
			else
			  p_part.push_back(ale->p_dip(g_id,gp_id,gpp_id));
			//cout << p_part.size() << " " ;
		      }
		  }
		}
	      else //at the root
		{
		  //root bipartition needs to be handled separately
		  map<set<long int>,int> bip_parts; // the map is here just for ordering the sets of clade ids. Each set only has 2 elements.
			for (map <long int,scalar_type> :: iterator it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
		    { //Going through all possible roots
				long int gp_id=(*it).first;
				boost::dynamic_bitset<> gamma=ale->id_sets.at(gp_id);
				boost::dynamic_bitset<> not_gamma = ~gamma;
                not_gamma[0] = 0;
				/*for (set<int>::iterator st=ale->Gamma.begin();st!=ale->Gamma.end();st++)
					if (gamma.count(*st)==0)
						not_gamma.insert(*st); //Building a function for constructing not_gamma would be useful!
                 */
				long int gpp_id = ale->set_ids.at(not_gamma);
				set <long int> parts;
				parts.insert(gp_id);
				parts.insert(gpp_id);
				bip_parts[parts]=1; //1 is a default value of no interest
				/*gamma.clear();
				not_gamma.clear();*/
		    }
		  for (map<set<long int>,int> :: iterator kt = bip_parts.begin();kt!=bip_parts.end();kt++)
		    {
		      vector <long int> parts;
		      for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
		      long int gp_id=parts[0];
		      long int gpp_id=parts[1];	    
		      gp_ids.push_back(gp_id);
		      gpp_ids.push_back(gpp_id);
		      if (ale->Bip_counts.at(gp_id)<=scalar_parameter.at("min_bip_count") and not ale->Gamma_size<4)
			p_part.push_back(0);	      
		      else
			p_part.push_back(ale->p_bip(gp_id));	      
		    }
		  bip_parts.clear();
		}
             //Now we have filled the vectors gp_ids, gpp_ids and p_part: we know all resolutions of clade g_id with the associated probability.
			
	      int N_parts=gp_ids.size(); //Number of resolutions of clade g_id.
	      if (!1) { //N_parts >= num_threads) { // It makes sense to do further parallelization but that slows things down!
		//iterate over all postions along S
		for (int rank=0;rank<last_rank;rank++)
		  {
		    int n=time_slices[rank].size();
		    for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
		      {
			//######################################################################################################################
			//#########################################INNNER LOOP##################################################################
			//######################################################################################################################
                         
			scalar_type t=time_slice_times[rank][t_i];
			scalar_type tpdt;//,tpdt_nl;
			int tpdt_rank, tpdt_t_i;
			if ( t_i <  (int)time_slice_times[rank].size()-1 ){
			  tpdt=time_slice_times[rank][t_i+1];
			  tpdt_rank = rank;
			  tpdt_t_i = t_i+1;
			}
			else if (rank<last_rank-1){
			  tpdt=time_slice_times[rank+1][0];
			  tpdt_rank = rank+1;
			  tpdt_t_i = 0;
			}
			else{
			  //top of root stem
			  tpdt=t_begin[time_slices[rank][0]]; //PBM PBM PBM
			  tpdt_rank = rank;  //PBM PBM PBM
			  tpdt_t_i = 0;	 //PBM PBM PBM		
			}
			
			bool tpdt_nl_is_t = true;
			if (scalar_parameter["event_node"]==1 and false)
			  ;//tpdt_nl=t;			
			else{
			  //tpdt_nl=tpdt;
			  tpdt_nl_is_t = false;
			}


		        //root
			scalar_type Delta_t=tpdt-t;
			//scalar_type N=vector_parameter["N"][rank];
			scalar_type Delta_bar=vector_parameter["Delta_bar"][rank];
			//scalar_type Lambda_bar=vector_parameter["Lambda_bar"][rank];
			//OMG
			//scalar_type p_Delta_bar=1-exp(-Delta_bar/N*Delta_t);			     
			scalar_type p_Delta_bar=Delta_bar*Delta_t;			     
			scalar_type Ebar=Ee[-1][t];
                         
			//boundaries for branch alpha virtual branch  
			//boundary at present
			if (t==0) {
			  qvec[g_id+1][rank][t_i][alpha]=0;
			  //q[g_id][t][alpha]=0;
			}
			//boundary between slice rank and rank-1 slice is trivial	
			;//q[g_id][t][alpha]=q[g_id][t][alpha];	  
                         //boundaries for branch alpha virtual branch.  
			if(1)
			  {
#pragma omp parallel for schedule(dynamic,1) //p9
			    for (int branch_i=0;branch_i<n;branch_i++)
			      {	  
				{
				  int e = time_slices[rank][branch_i];
                                     
				  //boundaries for branch e
				  //boundary at present
				  if (t==0)
				    {
				      if (is_a_leaf && extant_species[e]==gid_sps[g_id]) {	
					qvec[g_id+1][rank][t_i][e]=1;	
					//q[g_id][t][e]=1;
				      }
				      else {
					qvec[g_id+1][rank][t_i][e]=0;
					//q[g_id][t][e]=0;
				      }
				    }
				  //boundary between slice rank and rank-1
				  else if (t_i==0)
				    {
				      //terminating branch is last in time_slices and defines a represented speciation 
				      if (branch_i==n-1 && rank>0)
					{
					  int f=daughters[e][0];
					  int g=daughters[e][1];
					  scalar_type Eft=Ee[f][t];
					  scalar_type Egt=Ee[g][t];
                                             
					  scalar_type q_sum=0;
					  //q[g_id][t][e]=0;
                                             
					  scalar_type SL_fLg=qvec[g_id+1][rank][t_i][f]*Egt;
					  scalar_type SL_Lfg=qvec[g_id+1][rank][t_i][g]*Eft;
					  //scalar_type SL_fLg=q[g_id][t][f]*Egt;
					  //scalar_type SL_Lfg=q[g_id][t][g]*Eft;
					  //SL EVENT
					  q_sum+=SL_fLg+SL_Lfg;
					  //q[g_id][t][e]=q[g_id][t][f]*Egt + q[g_id][t][g]*Eft;
					  //SL.
                                             
					  //non-leaf directed partition
					  if (not is_a_leaf)
					    for (int i=0;i<N_parts;i++)
					      {	
						long int gp_id=gp_ids[i];
						long int gpp_id=gpp_ids[i];	    
						scalar_type pp=p_part[i];
						scalar_type S_pf_ppg=qvec[gp_id+1][rank][t_i][f]*qvec[gpp_id+1][rank][t_i][g]*pp;
						scalar_type S_ppf_pg=qvec[gpp_id+1][rank][t_i][f]*qvec[gp_id+1][rank][t_i][g]*pp;
						//scalar_type S_pf_ppg=q[gp_id][t][f]*q[gpp_id][t][g]*pp;
						//scalar_type S_ppf_pg=q[gpp_id][t][f]*q[gp_id][t][g]*pp;
						//S EVENT
						//q[g_id][t][e]+=q[gp_id][t][f]*q[gpp_id][t][g] +q[gpp_id][t][f]*q[gp_id][t][g];
						q_sum+= S_pf_ppg + S_ppf_pg;
						//S.
                                                     
					      }
					  qvec[g_id+1][rank][t_i][e]=q_sum;
					  //q[g_id][t][e]=q_sum; 
                                             
					}
				      //branches that cross to next time slice  
				      else
					{
					  //trivial
					  ;//q[g_id][t][e]=q[g_id][t][e];
					}			  
				    }		   
				  //boundaries for branch e.
				}
			      }
			  }
                         
			if(1)
			  {
			    //events within slice rank at time t on alpha virtual branch
			    scalar_type G_bar=Ge[-1][t];//exp(-(Delta_bar*(n-N)/N+Lambda_bar)*Delta_t );
			    qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]=0;
			    //q[g_id][tpdt][alpha]=0;
			    scalar_type q_sum=0;
			    scalar_type q_sum_nl=0;
#pragma omp parallel for schedule(dynamic,1)  reduction(+:q_sum_nl) //p10
			    for (int branch_i=0;branch_i<n;branch_i++)			  
			      {                
				int e = time_slices[rank][branch_i];		
				scalar_type tau_e=vector_parameter["tau"][e];
				//G_bar*=exp(- tau_e*Delta_t);
                                 
				//scalar_type p_Ntau_e=1-exp(-N*tau_e*Delta_t);
				//OMG
				scalar_type p_Ntau_e=tau_e*Delta_t;
				//non-leaf directed partition
				if (not is_a_leaf)
				  for (int i=0;i<N_parts;i++)
				    {	
				      long int gp_id=gp_ids[i];
				      long int gpp_id=gpp_ids[i];	    
				      scalar_type pp=p_part[i];				    
				      scalar_type T_ep_app=p_Ntau_e*qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]*pp;
				      scalar_type T_ap_epp=p_Ntau_e*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]*pp;
				      //scalar_type T_ep_app=p_Ntau_e*q[gp_id][t][e]*q[gpp_id][t][alpha]*pp;
				      //scalar_type T_ap_epp=p_Ntau_e*q[gp_id][t][alpha]*q[gpp_id][t][e]*pp;
				      //T EVENT
				      q_sum_nl+=T_ep_app+T_ap_epp;
				      //q[g_id][tpdt][alpha]+=p_Ntau_e*(q[gp_id][t][e]*q[gpp_id][t][alpha]+q[gp_id][t][alpha]*q[gpp_id][t][e]);
				      //T.
                                         
				    }
                                 
			      }
			    //non-leaf directed partition
			    if (not is_a_leaf)
			      for (int i=0;i<N_parts;i++)
				{	
				  long int gp_id=gp_ids[i];
				  long int gpp_id=gpp_ids[i];	    
				  scalar_type pp=p_part[i];
				  scalar_type Sb=p_Delta_bar*(2*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha])*pp;
				  //scalar_type Sb=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha])*pp;
				  //S_bar EVENT
				  q_sum_nl+=Sb;
				  //q[g_id][tpdt][alpha]+=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha]);
				  //S_bar.
                                     
				}			    
			    if (tpdt_nl_is_t)
			      qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum_nl;
			    else
			      qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum_nl;			    
			    //q[g_id][tpdt_nl][alpha]+=q_sum_nl;
			    for (int branch_i=0;branch_i<n;branch_i++)			  
			      {
                                 
				int e = time_slices[rank][branch_i];		
				scalar_type tau_e=vector_parameter["tau"][e];
				//OMG
				//scalar_type p_Ntau_e=1-exp(-N*tau_e*Delta_t);
				scalar_type p_Ntau_e=tau_e*Delta_t;
				scalar_type TLb=p_Ntau_e*Ebar*qvec[g_id+1][rank][t_i][e];
				//scalar_type TLb=p_Ntau_e*Ebar*q[g_id][t][e];
				//TL_bar EVENT
				q_sum+=TLb;
				//q[g_id][tpdt][alpha]+=p_Ntau_e*Ebar*q[g_id][t][e];
				//TL_bar.
                                 
			      }
			    //0 EVENT
			    scalar_type empty=G_bar*qvec[g_id+1][rank][t_i][alpha]; 
			    //scalar_type empty=G_bar*q[g_id][t][alpha]; 
			    q_sum+=empty;
                             
			    //q[g_id][tpdt][alpha]+=G_bar*q[g_id][t][alpha];
			    //0.
			    //max
			    /*
                              if (max_term<empty) 
                              {
                              max_term=empty;
                              }
			    */
			    //max		    
			    qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum;	    
			    //q[g_id][tpdt][alpha]+=q_sum;
			    //events within slice rank at time t on alpha virtual branch.
			  }
			if(1)
			  {
			    //DOES PREVENT SEGFAULT
#pragma omp parallel for schedule(dynamic,1) //p1
			    for (int branch_i=0;branch_i<n;branch_i++)
			      {	
                                 
				int e = time_slices[rank][branch_i];
				scalar_type Get=Ge[e][t];
				scalar_type Eet=Ee[e][t];	
				scalar_type delta_e=vector_parameter["delta"][e];
				scalar_type p_delta_e=delta_e*Delta_t;
                                 
				//events within slice rank at time t on branch e 
#pragma omp critical 
				{                                     
				  qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=0;
				  //q[g_id][tpdt][e]=0;
				}
				scalar_type q_sum=0;
				scalar_type q_sum_nl=0;
                                 
				//non-leaf directed partition		   
                                 
				if (not is_a_leaf)
				  for (int i=0;i<N_parts;i++)
				    {	
				      long int gp_id=gp_ids[i];
				      long int gpp_id=gpp_ids[i];	    
				      scalar_type pp=p_part[i];
				      scalar_type qpe=qvec[gp_id+1][rank][t_i][e];
				      scalar_type qppe=qvec[gpp_id+1][rank][t_i][e];
				      //scalar_type qpe=q[gp_id][t][e];
				      //scalar_type qppe=q[gpp_id][t][e];
				      scalar_type Sb_pa_ppe= p_Delta_bar*qvec[gp_id+1][rank][t_i][alpha]*qppe*pp;
				      scalar_type Sb_pe_ppa= p_Delta_bar*qpe*qvec[gpp_id+1][rank][t_i][alpha]*pp;
				      //scalar_type Sb_pa_ppe= p_Delta_bar*q[gp_id][t][alpha]*qppe*pp;
				      //scalar_type Sb_pe_ppa= p_Delta_bar*qpe*q[gpp_id][t][alpha]*pp;
				      //S_bar EVENT
				      q_sum_nl+= Sb_pa_ppe + Sb_pe_ppa;
                                         
				      //q[g_id][tpdt][e]+=p_Delta_bar*(q[gp_id][t][alpha]*q[gpp_id][t][e]+q[gp_id][t][e]*q[gpp_id][t][alpha]);			  
				      //S_bar.
                                         
				      scalar_type D=p_delta_e*qpe*qppe*pp;
				      //D EVENT
				      q_sum_nl+= D;
                                         
				      //q[g_id][tpdt][e]+=p_delta_e*q[gp_id][t][e]*q[gpp_id][t][e];
				      //D.
                                         
				    }
                                 
				scalar_type SLb=p_Delta_bar*Eet*qvec[g_id+1][rank][t_i][alpha];
				//scalar_type SLb=p_Delta_bar*Eet*q[g_id][t][alpha];
				//SL_bar EVENT
				q_sum_nl+=SLb;
                                 
				//q[g_id][tpdt][e]+=p_Delta_bar*Eet*q[g_id][t][alpha];
				//SL_bar.
#pragma omp critical 
				{ 
				  if (tpdt_nl_is_t)
				    qvec[g_id+1][rank][t_i][e]+=q_sum_nl;
				  else
				    qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum_nl;				  
				  //q[g_id][tpdt_nl][e]+=q_sum_nl;
				  
				}
				scalar_type empty=Get*qvec[g_id+1][rank][t_i][e];
				//scalar_type empty=Get*q[g_id][t][e];
				//0 EVENT
				q_sum+=empty;
                                 
				//q[g_id][tpdt][e]=Get*q[g_id][t][e];
				//0.
#pragma omp critical 
				{
				  qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum;
				  //q[g_id][tpdt][e]+=q_sum;
				}
				//events within slice rank at time t on branch e.
                                 
			      }
			  }
			//######################################################################################################################
			//#########################################INNNER LOOP##################################################################
			//######################################################################################################################		    
		      }
		  }
                 
	      }
	      else { // It does not make sense to do further parallelization
		//iterate over all time slices along S, from the leaves to the top
		for (int rank=0;rank<last_rank;rank++)
		  {
		    int n=time_slices[rank].size(); //Number of branches going through slice
		    for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
		      { 
				  //Going through all subslices, and for each subslice considering the set of possible SDTL events.
			//######################################################################################################################
			//#########################################INNNER LOOP##################################################################
			//######################################################################################################################
                         
			scalar_type t=time_slice_times[rank][t_i]; //End time of the current subslice.
			scalar_type tpdt;//,tpdt_nl; //tpdt: beginning time of the current subslice. tpdt_nl=tpdt is for the event node.
			int tpdt_rank, tpdt_t_i;
			if ( t_i <  (int)time_slice_times[rank].size()-1 ){
			  tpdt=time_slice_times[rank][t_i+1];
			  tpdt_rank = rank;
			  tpdt_t_i = t_i+1;
			}
			else if (rank<last_rank-1){
			  tpdt=time_slice_times[rank+1][0];
			  tpdt_rank = rank+1;
			  tpdt_t_i = 0;
			}
			else{
			  //top of root stem
			  tpdt=t_begin[time_slices[rank][0]]; //PBM PBM PBM
			  tpdt_rank = rank;  //PBM PBM PBM
			  tpdt_t_i = 0;	 //PBM PBM PBM		
			}
			
			bool tpdt_nl_is_t = true;
			if (scalar_parameter["event_node"]==1 and false)
			  ;//tpdt_nl=t;			
			else{
			  //tpdt_nl=tpdt;
			  tpdt_nl_is_t = false;
			}
			 
			//root
			scalar_type Delta_t=tpdt-t; //size of the subslice
			//scalar_type N=vector_parameter["N"][rank];
			scalar_type Delta_bar=vector_parameter["Delta_bar"][rank]; //Parameter for the speciation rate in the Moran process
			//scalar_type Lambda_bar=vector_parameter["Lambda_bar"][rank];
			//OMG
			//scalar_type p_Delta_bar=1-exp(-Delta_bar/N*Delta_t);			     
			scalar_type p_Delta_bar=Delta_bar*Delta_t;			     
			scalar_type Ebar=Ee[-1][t];
                         
			//boundaries for branch alpha virtual branch  
			//boundary at present
			if (t==0) {
			  //#pragma omp critical 
			  {
			    qvec[g_id+1][rank][t_i][alpha]=0;
			    //q[g_id][t][alpha]=0;
			  }
			}
			//boundary between slice rank and rank-1 slice is trivial	
			;//q[g_id][t][alpha]=q[g_id][t][alpha];	  
                         //boundaries for branch alpha virtual branch.  
			if(1) //Only dealing with S and SL events
			  {
			    for (int branch_i=0;branch_i<n;branch_i++) //For all branches going through the current slice
			      {	  
				{
				  int e = time_slices[rank][branch_i]; //index of a particular edge
                                     
				  //boundaries for branch e
				  //boundary at present
				  if (t==0)
				    {
				      if (is_a_leaf && extant_species[e]==gid_sps[g_id]) {	
					qvec[g_id+1][rank][t_i][e]=1;	
					//q[g_id][t][e]=1;
				      }
				      else {
					qvec[g_id+1][rank][t_i][e]=0;
					//q[g_id][t][e]=0;
				      }
				    }
				  //boundary between slice rank and rank-1
				  else if (t_i==0)
				    {
				      //terminating branch is last in time_slices and defines a represented speciation 
				      if (branch_i==n-1 && rank>0)
					{
					  int f=daughters[e][0];
					  int g=daughters[e][1];
					  scalar_type Eft=Ee[f][t];
					  scalar_type Egt=Ee[g][t];
                                             
					  scalar_type q_sum=0;
					  //q[g_id][t][e]=0;
                                             
					  scalar_type SL_fLg=qvec[g_id+1][rank][t_i][f]*Egt;
					  scalar_type SL_Lfg=qvec[g_id+1][rank][t_i][g]*Eft;
					  //scalar_type SL_fLg=q[g_id][t][f]*Egt;
					  //scalar_type SL_Lfg=q[g_id][t][g]*Eft;
					  //SL EVENT
					  q_sum+=SL_fLg+SL_Lfg;
					  //q[g_id][t][e]=q[g_id][t][f]*Egt + q[g_id][t][g]*Eft;
					  //SL.
                                             
					  //non-leaf directed partition
					  if (not is_a_leaf)
					    for (int i=0;i<N_parts;i++)
					      {	
						long int gp_id=gp_ids[i];
						long int gpp_id=gpp_ids[i];	    
						scalar_type pp=p_part[i];
						scalar_type S_pf_ppg=qvec[gp_id+1][rank][t_i][f]*qvec[gpp_id+1][rank][t_i][g]*pp;
						scalar_type S_ppf_pg=qvec[gpp_id+1][rank][t_i][f]*qvec[gp_id+1][rank][t_i][g]*pp;
						//scalar_type S_pf_ppg=q[gp_id][t][f]*q[gpp_id][t][g]*pp;
						//scalar_type S_ppf_pg=q[gpp_id][t][f]*q[gp_id][t][g]*pp;
						//S EVENT
						//q[g_id][t][e]+=q[gp_id][t][f]*q[gpp_id][t][g] +q[gpp_id][t][f]*q[gp_id][t][g];
						q_sum+= S_pf_ppg + S_ppf_pg;
						//S.
                                                     
					      }
					  qvec[g_id+1][rank][t_i][e]=q_sum; 
					  //q[g_id][t][e]=q_sum; 
                                             
					}
				      //branches that cross to next time slice  
				      else
					{
					  //trivial
					  ;//q[g_id][t][e]=q[g_id][t][e];
					}			  
				    }		   
				  //boundaries for branch e.
				}
			      }
			  }
                         
			if(1) //Only dealing with T events, S_bar, TL_bar and naught events?
			  {
			    //events within slice rank at time t on alpha virtual branch
			    scalar_type G_bar=Ge[-1][t];//exp(-(Delta_bar*(n-N)/N+Lambda_bar)*Delta_t );	
			    //#pragma omp critical 
			    {                                 
			      qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]=0;
			      //q[g_id][tpdt][alpha]=0;
			    }
			    scalar_type q_sum=0;
			    scalar_type q_sum_nl=0;
			    for (int branch_i=0;branch_i<n;branch_i++)			  
			      {
                                 
				int e = time_slices[rank][branch_i];		
				scalar_type tau_e=vector_parameter["tau"][e];
				//G_bar*=exp(- tau_e*Delta_t);
                                 
				//scalar_type p_Ntau_e=1-exp(-N*tau_e*Delta_t);
				//OMG
				scalar_type p_Ntau_e=tau_e*Delta_t;
				//non-leaf directed partition
				if (not is_a_leaf)
				  for (int i=0;i<N_parts;i++) //For each partition
				    {	
				      long int gp_id=gp_ids[i];
				      long int gpp_id=gpp_ids[i];	    
				      scalar_type pp=p_part[i];				    
				      scalar_type T_ep_app=p_Ntau_e*qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]*pp;
				      scalar_type T_ap_epp=p_Ntau_e*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]*pp;		    
				      //scalar_type T_ep_app=p_Ntau_e*q[gp_id][t][e]*q[gpp_id][t][alpha]*pp;
				      //scalar_type T_ap_epp=p_Ntau_e*q[gp_id][t][alpha]*q[gpp_id][t][e]*pp;
				      //T EVENT
				      q_sum_nl+=T_ep_app+T_ap_epp;
				      //q[g_id][tpdt][alpha]+=p_Ntau_e*(q[gp_id][t][e]*q[gpp_id][t][alpha]+q[gp_id][t][alpha]*q[gpp_id][t][e]);
				      //T.
                                         
				    }
                                 
			      }
			    //non-leaf directed partition
			    if (not is_a_leaf)
			      for (int i=0;i<N_parts;i++)
				{	
				  long int gp_id=gp_ids[i];
				  long int gpp_id=gpp_ids[i];	    
				  scalar_type pp=p_part[i];
				  scalar_type Sb=p_Delta_bar*(2*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha])*pp;
				  //scalar_type Sb=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha])*pp;
				  //S_bar EVENT Speciation in an unseen lineage?
				  q_sum_nl+=Sb;
				  //q[g_id][tpdt][alpha]+=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha]);
				  //S_bar.
                                     
				}	 
			    if (tpdt_nl_is_t)			      
			      qvec[g_id+1][rank][t_i][alpha]+=q_sum_nl;
			    else		      
			      qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum_nl;			      			    
			    //q[g_id][tpdt_nl][alpha]+=q_sum_nl;
			    for (int branch_i=0;branch_i<n;branch_i++)			  
			      {
                                 
				int e = time_slices[rank][branch_i];		
				scalar_type tau_e=vector_parameter["tau"][e];
				//OMG
				//scalar_type p_Ntau_e=1-exp(-N*tau_e*Delta_t);
				scalar_type p_Ntau_e=tau_e*Delta_t;
				scalar_type TLb=p_Ntau_e*Ebar*qvec[g_id+1][rank][t_i][e];
				//scalar_type TLb=p_Ntau_e*Ebar*q[g_id][t][e];
				//TL_bar EVENT
				q_sum+=TLb;
				//q[g_id][tpdt][alpha]+=p_Ntau_e*Ebar*q[g_id][t][e];
				//TL_bar.
                                 
			      }
			    //0 EVENT
			    scalar_type empty=G_bar*qvec[g_id+1][rank][t_i][alpha]; 
			    //scalar_type empty=G_bar*q[g_id][t][alpha]; 
			    q_sum+=empty;
                             
			    //q[g_id][tpdt][alpha]+=G_bar*q[g_id][t][alpha];
			    //0.
			    //max
			    /*
                              if (max_term<empty) 
                              {
                              max_term=empty;
                              }
			    */
			    //max		    
			    //#pragma omp critical 
			    {
                              qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum;
			      //q[g_id][tpdt][alpha]+=q_sum;
			    }
			    //events within slice rank at time t on alpha virtual branch.
			  }
				  
			if(1) //Only dealing with S_bar, D, naught events
			  {
			    for (int branch_i=0;branch_i<n;branch_i++)
			      {	
                                 
				int e = time_slices[rank][branch_i];
				scalar_type Get=Ge[e][t];
				scalar_type Eet=Ee[e][t];	
				scalar_type delta_e=vector_parameter["delta"][e];
				scalar_type p_delta_e=delta_e*Delta_t;
                                 
				//events within slice rank at time t on branch e 
				qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=0;
				//q[g_id][tpdt][e]=0;
				scalar_type q_sum=0;
				scalar_type q_sum_nl=0;
                                 
				//non-leaf directed partition		   
                                 
				if (not is_a_leaf)
				  for (int i=0;i<N_parts;i++)
				    {	
				      long int gp_id=gp_ids[i];
				      long int gpp_id=gpp_ids[i];	    
				      scalar_type pp=p_part[i];
				      scalar_type qpe=qvec[gp_id+1][rank][t_i][e];
				      scalar_type qppe=qvec[gpp_id+1][rank][t_i][e];
				      //scalar_type qpe=q[gp_id][t][e];
				      //scalar_type qppe=q[gpp_id][t][e];
				      scalar_type Sb_pa_ppe= p_Delta_bar*qvec[gp_id+1][rank][t_i][alpha]*qppe*pp;
				      scalar_type Sb_pe_ppa= p_Delta_bar*qpe*qvec[gpp_id+1][rank][t_i][alpha]*pp;
				      //scalar_type Sb_pa_ppe= p_Delta_bar*q[gp_id][t][alpha]*qppe*pp;
				      //scalar_type Sb_pe_ppa= p_Delta_bar*qpe*q[gpp_id][t][alpha]*pp;
				      //S_bar EVENT
				      q_sum_nl+= Sb_pa_ppe + Sb_pe_ppa;
                                         
				      //q[g_id][tpdt][e]+=p_Delta_bar*(q[gp_id][t][alpha]*q[gpp_id][t][e]+q[gp_id][t][e]*q[gpp_id][t][alpha]);			  
				      //S_bar.
                                         
				      scalar_type D=p_delta_e*qpe*qppe*pp;
				      //D EVENT
				      q_sum_nl+= D;
                                         
				      //q[g_id][tpdt][e]+=p_delta_e*q[gp_id][t][e]*q[gpp_id][t][e];
				      //D.
                                         
				    }
                                 
				scalar_type SLb=p_Delta_bar*Eet*qvec[g_id+1][rank][t_i][alpha];                                 
				//scalar_type SLb=p_Delta_bar*Eet*q[g_id][t][alpha];
				//SL_bar EVENT
				q_sum_nl+=SLb;
                                 
				//q[g_id][tpdt][e]+=p_Delta_bar*Eet*q[g_id][t][alpha];
				//SL_bar.  
				if (tpdt_nl_is_t)
				  qvec[g_id+1][rank][t_i][e]+=q_sum_nl;
				else
				  qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum_nl;				
				//q[g_id][tpdt_nl][e]+=q_sum_nl;
				scalar_type empty=Get*qvec[g_id+1][rank][t_i][e];
				//scalar_type empty=Get*q[g_id][t][e];
				//0 EVENT
				q_sum+=empty;
                                 
				//q[g_id][tpdt][e]=Get*q[g_id][t][e];
				//0.
				qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum;
				//q[g_id][tpdt][e]+=q_sum;
				//events within slice rank at time t on branch e.
                                 
			      }
			  }
			//######################################################################################################################
			//#########################################INNNER LOOP##################################################################
			//######################################################################################################################		    
		      }
		  }
	      }
             
	      gp_ids.clear();
	      gpp_ids.clear();
	      p_part.clear();
	      //scalar_type tnow=omp_get_wtime();//t->elapsed();
	      //n cout <<  "Outer loop parallelization: Nb Partitions: "<<N_parts << " iteration duration: "<<(tnow-tatom) << endl; ;
	      //tatom=tnow;
	    }
         
	}     
     
      }

      //scalar_type tnow=omp_get_wtime();//t->elapsed();
      //n cout << endl << it2->first  << " "<< siz << " "<<tnow-told << " " <<(tnow-told)/siz << endl;
      //told=tnow;

    }
  //n cout << "Entire loop duration: "<< omp_get_wtime() -toldest<< endl;
  //cout << "end loop" << endl;
  scalar_type root_norm=0;
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{	  
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {	      
	      root_norm+=1;	      
	    }	     
	  root_norm+=1;	      
	}
    }

  scalar_type root_sum=0;
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{
	  //scalar_type t=time_slice_times[rank][t_i];		
	  
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {
	      int e = time_slices[rank][branch_i];
	      root_sum+=qvec[0][rank][t_i][e]/root_norm;	   
	      //root_sum+=q[-1][t][e]/root_norm;	      
	    }	     
	  root_sum+=qvec[0][rank][t_i][alpha]/root_norm;	 
	  //root_sum+=q[-1][t][alpha]/root_norm;	      
	}
    }
  /*test
    for (int rank=0;rank<last_rank;rank++)
    {
    int n=time_slices[rank].size();
    for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
    {
    scalar_type t=time_slice_times[rank][t_i];
		
    cout << rank << " " << n << " " << t << " " << q[tmp_g_id][t][alpha]<< endl;
	  
    for (int branch_i=0;branch_i<n;branch_i++)
    {	    
    int branch = time_slices[rank][branch_i];
    scalar_type tmp_q=q[tmp_g_id][t][branch];
    Node * tmp_node = id_nodes[branch];
    //cout << branch << " " << t_i << " " << rank <<" "<< tmp_node << endl;
    stringstream out;
    string name = (* (dynamic_cast<const BppString *>(tmp_node->getBranchProperty("ID")))).toSTL();
    if (tmp_q>0)
    out << log(tmp_q) ;
    //out << Ge[branch][t];
    tmp_node->setBranchProperty("ID",BppString(name+out.str().substr(0,4)+"|"));	      
    //if (tmp_node->isLeaf())
    //tmp_node->setName(tmp_node->getName()+out.str().substr(0,4)+"|");
    }
    }
    }
    cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID") << endl;
    for (map <Node *,int >::iterator it=node_ids.begin();it!=node_ids.end();it++ )
    (*it).first->setBranchProperty("ID",BppString(""));
  */
  //test

  //del-locs
  g_ids.clear();
  g_id_sizes.clear();

  return root_sum;	
}

void exODT_model::calculate_EGb()
{

  for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ee.begin();it!=Ee.end();it++)//del_loc
    (*it).second.clear();
  Ee.clear();
  for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ge.begin();it!=Ge.end();it++)//del_loc
    (*it).second.clear();
  Ge.clear();


  map<int,scalar_type> Ee_y;//del-loc
  map<int,scalar_type> Ge_y;//del-loc
  map<int,scalar_type> E_k1,E_k2,E_k3,E_k4;//del-loc
  map<int,scalar_type> G_k1,G_k2,G_k3,G_k4;//del-loc

  for (int rank=0;rank<last_rank;rank++) 
    for (int tsi=0;tsi<(int)time_slice_times[rank].size();tsi++)
      {
	map<int,map <scalar_type,scalar_type> > y_E,y_G;//del-loc
	map<int,map <int,scalar_type> > iy_E,iy_G;//del-loc

	scalar_type t_b;
	if (tsi==(int)time_slice_times[rank].size()-1)
	  t_b = time_slice_begins[rank];
	else
	  t_b = time_slice_times[rank][tsi+1];
	scalar_type t_e;
	if (tsi==0)
	  {
	    if (rank>0 )
	      t_e = time_slice_begins[rank-1];
	    else
	      t_e = 0;
	  }
	else
	  {
	    t_e=time_slice_times[rank][tsi];
	  }
	scalar_type N=vector_parameter["N"][rank];

	scalar_type ni=time_slices[rank].size();
	scalar_type Delta_bar=vector_parameter["Delta_bar"][rank];//1
	scalar_type Lambda_bar=vector_parameter["Lambda_bar"][rank]*N/(N-ni);;
      
	scalar_type t=t_e;
	scalar_type tpdt=t_b;
	scalar_type h=(tpdt-t)/scalar_parameter["DD"];
	//scalar_type ti=t;

	scalar_type h_lambda_avg=h*scalar_parameter["lambda_avg"]; 
	scalar_type h_delta_avg=h*scalar_parameter["delta_avg"];	      
	scalar_type h_tau_avg=h*scalar_parameter["tau_avg"]*(N-ni)/(N-1)*N;	      
	scalar_type h_Delta_bar=h*Delta_bar;
	scalar_type h_Lambda_bar=h*Lambda_bar;


	for (int ii=0;ii<scalar_parameter["DD"];ii++)
	  {

	    //intial conditions
	    if (ii==0)
	      {
		if ( t==0)
		  Ee[-1][t]=1;
		//trivial else Ee[-1][t]=Ee[-1][t];
	      
		//y_E[-1][t]=Ee[-1][t];
		iy_E[-1][ii]=Ee[-1][t];

		//Ee_y[-1]=y_E[-1][t];
		Ee_y[-1]=iy_E[-1][ii];
	      
		Ge_y[-1]=1;

		//y_G[-1][t]=1;	   
		iy_G[-1][ii]=1;	   

	      }
	  
	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		if (ii==0)
		  {
		    if ( t==0)
		      { 
			Ee[e][t]=0;
		      }
		    else if (t==t_end[e])
		      {
			int f=daughters[e][0];int g=daughters[e][1];
			Ee[e][t]=Ee[f][t]*Ee[g][t];	   
		      }
		    //trivial else{Ee[e][t]=Ee[e][t];}
		    //y_E[e][t]=Ee[e][t];
		    iy_E[e][ii]=Ee[e][t];

		    //Ee_y[e]=y_E[e][t];
		    Ee_y[e]=iy_E[e][ii];

		    Ge_y[e]=1;

		    //y_G[e][t]=1;
		    iy_G[e][ii]=1;

		  }
	      }
	    // RK4: 4th order Runge-Kutta for y'=f(y) 
	    // k1 = f(y[n])	      
	    E_k1[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k1[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];	       
		E_k1[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k1[-1]-=h_tau_f*(1-Ee_y[f])* Ge_y[-1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda; 
		scalar_type h_delta=h*delta;	      

		// k1 = f(y[n])	      
		E_k1[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k1[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];

	      }
	    // k2 = f(y[n]+h/2 k1)

	    //Ee_y[-1]=y_E[-1][ti]+1/2.* E_k1[-1];	      
	    Ee_y[-1]=iy_E[-1][ii]+1/2.* E_k1[-1];	      
	    //Ge_y[-1]=y_G[-1][ti]+1/2.* G_k1[-1];	      
	    Ge_y[-1]=iy_G[-1][ii]+1/2.* G_k1[-1];	      

	    E_k2[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k2[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];	       
		E_k2[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k2[-1]-=h_tau_f*(1-Ee_y[f])* Ge_y[-1];
	      }
	  	  
	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda; 
		scalar_type h_delta=h*delta;	      
		
		// k2 = f(y[n]+h/2 k1)
		//Ee_y[e] =y_E[e][ti]+1/2. * E_k1[e];
		Ee_y[e] =iy_E[e][ii]+1/2. * E_k1[e];
		//Ge_y[e] =y_G[e][ti]+1/2. * G_k1[e];
		Ge_y[e] =iy_G[e][ii]+1/2. * G_k1[e];
		
		E_k2[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k2[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];
	      }

	    // k3 = f(y[n]+h/2 k2)
	    //Ee_y[-1]=y_E[-1][ti]+1/2.* E_k2[-1];
	    Ee_y[-1]=iy_E[-1][ii]+1/2.* E_k2[-1];
	    //Ge_y[-1]=y_G[-1][ti]+1/2.* G_k2[-1];
	    Ge_y[-1]=iy_G[-1][ii]+1/2.* G_k2[-1];

	    E_k3[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k3[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];	       
		E_k3[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k3[-1]-=h_tau_f*(1-Ee_y[f])* Ge_y[-1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda; 
		scalar_type h_delta=h*delta;	      
 
		// k3 = f(y[n]+h/2 k2)
		//Ee_y[e] =y_E[e][ti]+1/2. * E_k2[e];
		Ee_y[e] =iy_E[e][ii]+1/2. * E_k2[e];
		//Ge_y[e] =y_G[e][ti]+1/2. * G_k2[e];
		Ge_y[e] =iy_G[e][ii]+1/2. * G_k2[e];

		E_k3[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k3[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];
	      }

	    // k4 = f(y[n]+h k3)
	    //Ee_y[-1]=y_E[-1][ti]+1* E_k3[-1];
	    Ee_y[-1]=iy_E[-1][ii]+1* E_k3[-1];
	    //Ge_y[-1]=y_G[-1][ti]+1* G_k3[-1];
	    Ge_y[-1]=iy_G[-1][ii]+1* G_k3[-1];

	    E_k4[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k4[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];	       
		E_k4[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k4[-1]-=h_tau_f*(1-Ee_y[f]) *Ge_y[-1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda; 
		scalar_type h_delta=h*delta;	      

		// k4 = f(y[n]+h k3)
		//Ee_y[e] =y_E[e][ti]+1 * E_k3[e];
		Ee_y[e] =iy_E[e][ii]+1 * E_k3[e];

		//Ge_y[e] =y_G[e][ti]+1 * G_k3[e];
		Ge_y[e] =iy_G[e][ii]+1 * G_k3[e];

		E_k4[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k4[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e])  + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];
	      }	  
	    // y[n+1] = y[n] + h/6 (k1 + 2 k2 + 2 k3 + k4) 
	    //y_E[-1][ti+h]=Ee_y[-1] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);
	    iy_E[-1][ii+1]=Ee_y[-1] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);	    
	    
	    ///*
	    if (ii==0)
	      iy_E[-1][ii+1]=Ee[-1][t] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);
	    else
	      iy_E[-1][ii+1]=iy_E[-1][ii] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);
	    //*/
	    //y_G[-1][ti+h]=Ge_y[-1] + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);
	    //iy_G[-1][ii+1]=Ge_y[-1] + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);

	    
	    if (ii==0)
	      iy_G[-1][ii+1]=1 + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);
	    else
	      iy_G[-1][ii+1]=iy_G[-1][ii] + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);
	    

	    if (ii==scalar_parameter["DD"]-1)
	      {
		//Ee[-1][tpdt]=y_E[-1][ti+h];
		Ee[-1][tpdt]=iy_E[-1][ii+1];

		//Ge[-1][t]=y_G[-1][ti+h];
		Ge[-1][t]=iy_G[-1][ii+1];

	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		// y[n+1] = y[n] + h/6 (k1 + 2 k2 + 2 k3 + k4) 
		//y_E[e][ti+h]=Ee_y[e] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);
		iy_E[e][ii+1]=Ee_y[e] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);		

		///*
		if (ii==0)
		  iy_E[e][ii+1]=Ee[e][t] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);
		else
		  iy_E[e][ii+1]=iy_E[e][ii] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);
		//*/

		//y_G[e][ti+h]=Ge_y[e] + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);	      		
		iy_G[e][ii+1]=Ge_y[e] + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);	      

		if (ii==0)
		  iy_G[e][ii+1]=1 + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);	      
		else		  
		  iy_G[e][ii+1]=iy_G[e][ii] + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);	      

		if (ii==scalar_parameter["DD"]-1)
		  {
		    //Ee[e][tpdt]=y_E[e][ti+h];
		    Ee[e][tpdt]=iy_E[e][ii+1];
		    //Ge[e][t]=y_G[e][ti+h];
		    Ge[e][t]=iy_G[e][ii+1];
		  }
	      }
	    //ti=ti+h;	  

	  }      
      }
  //del-locs
  Ee_y.clear();
  Ge_y.clear();
  E_k1.clear();E_k2.clear();E_k3.clear();E_k4.clear();
  G_k1.clear();G_k2.clear();G_k3.clear();G_k4.clear();
  /*  
  for (map<int,map <scalar_type,scalar_type> >::iterator it=y_E.begin();it!=y_E.end();it++)
    (*it).second.clear();
  y_E.clear();
  for (map<int,map <scalar_type,scalar_type> >::iterator it=y_G.begin();it!=y_G.end();it++)
    (*it).second.clear();
  y_G.clear();	 	   
  */
}
