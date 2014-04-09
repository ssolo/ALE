#include "exODT.h"
using namespace std;
using namespace bpp;

//The current lowmem=true method uses sample(true) cf. sample.cpp. 
//The general structure of the calculation, and lot of the code, is the same as p(ale) cf. model.cpp.
pair<string,scalar_type> exODT_model::p_MLRec(approx_posterior *ale, bool lowmem)
{

  ale_pointer=ale;
  //directed partitions and thier sizes
  vector <long int>  g_ids;//del-loc
  vector <long int>  g_id_sizes;//del-loc
  for (map <int, vector <long int > > :: iterator it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	g_ids.push_back((*jt));
	g_id_sizes.push_back((*it).first);
      }
  //root biprartitino needs to be handled seperatly
  g_ids.push_back(-1);
  g_id_sizes.push_back(ale->Gamma_size);

  // gene<->species mapping
  //vector<vector<vector<map<int, scalar_type> > > > qvec;
  qvec.clear();//hope this doesn't leak..
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
		    
		    //qvec[g_id+1][rank][t_i][e]=0; //initializing q with 0s for each clade of the ale, each time subslice, each branch
		    vbranch[e]=0;
		  }
		//qvec[g_id+1][rank][t_i][alpha]=0; //I don't understand the purpose of that: alpha=-1, so it's already 0, right? 
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

  //p_parts is filled up with CCPs
  for (int i=0;i<(int)g_ids.size();i++)
    {	
      // directed partition (dip) gamma's id  
      bool is_a_leaf=false;
      long int g_id=g_ids[i];	
      if (g_id_sizes[i]==1)
	is_a_leaf=true;
      
      vector <long int> gp_ids;//del-loc
      vector <long int> gpp_ids;//del-loc
      vector <scalar_type> p_part;//del-loc
      if (g_id!=-1)
	for (unordered_map< pair<long int, long int>,scalar_type> :: iterator kt = ale->Dip_counts[g_id].begin(); kt != ale->Dip_counts[g_id].end(); kt++)
	  {	  
	    pair<long int, long int> parts = (*kt).first;
	    long int gp_id=parts.first;
	    long int gpp_id=parts.second;
	    gp_ids.push_back(gp_id);
	    gpp_ids.push_back(gpp_id);
	    if (ale->Bip_counts[g_id]<=scalar_parameter["min_bip_count"])
	      p_part.push_back(0);
	    else
	      {
		p_part.push_back(ale->p_dip(g_id,gp_id,gpp_id));
	      }
	  }
      else
	{
	  //root biprartition needs to be handled seperatly
	  map<set<long int>,int> bip_parts;
	  for (map <long int,scalar_type> :: iterator it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
	    {
	      long int gp_id=(*it).first;
	      boost::dynamic_bitset<>  gamma = ale->id_sets[gp_id];
	      boost::dynamic_bitset<>  not_gamma = ~gamma;
	      not_gamma[0] = 0;
	      /*	      for (set<int>::iterator st=ale->Gamma.begin();st!=ale->Gamma.end();st++)
			      if (gamma.count(*st)==0)
			      not_gamma.insert(*st);*/
	      /*  for (auto i = 0; i < ale->nbint; ++i) {
		  not_gamma[i] = 0;
		  }
		  BipartitionTools::bitNot(not_gamma, gamma, ale->nbint);*/
	      long int gpp_id = ale->set_ids[not_gamma];
	      set <long int> parts;
	      parts.insert(gp_id);
	      parts.insert(gpp_id);
	      bip_parts[parts]=1;
	      /* gamma.clear();
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
	      if (ale->Bip_counts[gp_id]<=scalar_parameter["min_bip_count"] and not ale->Gamma_size<4)
		p_part.push_back(0);	      
	      else
		p_part.push_back(ale->p_bip(gp_id));	      
	    }
	  bip_parts.clear();
	}
      int N_parts=gp_ids.size();
  
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
	      scalar_type tpdt;
	      int tpdt_rank, tpdt_t_i;       			
	      if ( t_i < (int)time_slice_times[rank].size()-1 )
		{
		  tpdt=time_slice_times[rank][t_i+1];
		  tpdt_rank=rank;
		  tpdt_t_i=t_i+1;
		}
	      else if (rank<last_rank-1)
		{
		  tpdt=time_slice_times[rank+1][0];
		  tpdt_rank=rank+1;
		  tpdt_t_i=0;
		}
	      else
		//top of root ste
		{
		  tpdt=t_begin[time_slices[rank][0]];
		  tpdt_rank=rank;
		  tpdt_t_i=0;
		}


	      //root
	      scalar_type Delta_t=tpdt-t;
	      scalar_type Delta_bar=vector_parameter["Delta_bar"][rank];
	      scalar_type p_Delta_bar=Delta_bar*Delta_t;
			      
	      scalar_type Ebar=Ee[-1][t];;
	      if(1)
		{
		  for (int branch_i=0;branch_i<n;branch_i++)
		    {	    
		      int e = time_slices[rank][branch_i];

		      //boundaries for branch e
		      //boundary at present
		      if (t==0)
			{
			  if (is_a_leaf && extant_species[e]==gid_sps[g_id])			  
			    qvec[g_id+1][rank][t_i][e]=1;
			  else
			    qvec[g_id+1][rank][t_i][e]=0;
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
			    
			      //sum scalar_type q_sum=0;
			      //qvec[g_id+1][rank][t_i][e]=0;
			      //max
			      scalar_type max_term=0;
			      //max
			    
			      scalar_type SL_fLg=qvec[g_id+1][rank][t_i][f]*Egt;
			      scalar_type SL_Lfg=qvec[g_id+1][rank][t_i][g]*Eft;
			      //SL EVENT
			      // q_sum+=SL_fLg+SL_Lfg;
			      //qvec[g_id+1][rank][t_i][e]=qvec[g_id+1][rank][t_i][f]*Egt + qvec[g_id+1][rank][t_i][g]*Eft;
			      //SL.
			      //max
			      if (max_term<SL_fLg) 
				{
				  max_term= SL_fLg;
				}
			      if (max_term<SL_Lfg) 
				{
				  max_term= SL_Lfg;
				}
			      //max

			      //non-leaf directed partition
			      if (not is_a_leaf)
				for (int i=0;i<N_parts;i++)
				  {	
				    long int gp_id=gp_ids.at(i);
				    long int gpp_id=gpp_ids.at(i);	    
				    scalar_type pp=p_part.at(i);
				    scalar_type S_pf_ppg=qvec[gp_id+1][rank][t_i][f]*qvec[gpp_id+1][rank][t_i][g]*pp;
				    scalar_type S_ppf_pg=qvec[gpp_id+1][rank][t_i][f]*qvec[gp_id+1][rank][t_i][g]*pp;
				    //S EVENT
				    //qvec[g_id+1][rank][t_i][e]+=qvec[gp_id+1][rank][t_i][f]*qvec[gpp_id+1][rank][t_i][g] +qvec[gpp_id+1][rank][t_i][f]*qvec[gp_id+1][rank][t_i][g];
				    //sum q_sum+= S_pf_ppg + S_ppf_pg;
				    //S.
				    //max
				    if (max_term<S_pf_ppg) 
				      {
					max_term= S_pf_ppg;
				      }
				    if (max_term<S_ppf_pg) 
				      {
					max_term=S_ppf_pg;
				      }
				    //max

				  }

			      //sum qvec[g_id+1][rank][t_i][e]=q_sum; 
			      qvec[g_id+1][rank][t_i][e]=max_term; 
			    }
			  //branches that cross to next time slice  
			  else
			    {
			      //trivial
			      ;//qvec[g_id+1][rank][t_i][e]=qvec[g_id+1][rank][t_i][e];
			    }			  
			}		   
		      //boundaries for branch e.
		    }
		}	    
	      if(1)
		{
		  //boundaries for branch alpha virtual branch  
		  //boundary at present
		  if (t==0)
		    qvec[g_id+1][rank][t_i][alpha]=0;
		  //boundary between slice rank and rank-1 slice is trivial	
		  ;//qvec[g_id+1][rank][t_i][alpha]=qvec[g_id+1][rank][t_i][alpha];	  
		  //boundaries for branch alpha virtual branch.  

		  //events within slice rank at time t on alpha virtual branch
		  scalar_type G_bar=Ge[-1][t];
		  qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]=0;
		  //sum scalar_type q_sum=0;
		  //max
		  scalar_type max_term=0;
		  //max
		  for (int branch_i=0;branch_i<n;branch_i++)			  
		    {
		      int e = time_slices[rank][branch_i];		
		      scalar_type tau_e=vector_parameter["tau"][e];
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
			    //Tb EVENT
			    //sum q_sum+=T_ep_app+T_ap_epp;
			    //qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=p_Ntau_e*(qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]+qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]);
			    //Tb.
			    //max
			    if (max_term<T_ep_app) 
			      {
				max_term=T_ep_app;
			      }
			    if (max_term<T_ap_epp) 
			      {
				max_term=T_ap_epp;
			      }
			    //max

			  }
		    }
		  //non-leaf directed partition
		  if (not is_a_leaf)
		    for (int i=0;i<N_parts;i++)
		      {	
			long int gp_id=gp_ids.at(i);
			long int gpp_id=gpp_ids[i];	    
			scalar_type pp=p_part[i];
			scalar_type Sb=p_Delta_bar*(2*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha])*pp;
			//S_bar EVENT
			//sum q_sum+=Sb;
			//qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=p_Delta_bar*(2*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha]);
			//S_bar.
			//max
			if (max_term<Sb) 
			  {
			    max_term=Sb;
			  }
			//max

		      }	    
	      
		  for (int branch_i=0;branch_i<n;branch_i++)			  
		    {
		      int e = time_slices[rank][branch_i];		
		      scalar_type tau_e=vector_parameter["tau"][e];
		      scalar_type p_Ntau_e=tau_e*Delta_t;
		      scalar_type TLb=p_Ntau_e*Ebar*qvec[g_id+1][rank][t_i][e];
		      //TL_bar EVENT
		      //sum q_sum+=TLb;
		      //qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=p_Ntau_e*Ebar*qvec[g_id+1][rank][t_i][e];
		      //TL_bar.
		      //max
		      if (max_term<TLb) 
			{
			  max_term=TLb;
			}
		      //max		    
		    }
	      
		  //sum qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum_nl;
		  if (qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]<max_term)		      
		    {
		      qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]=max_term;
		    }
	      
		  //0 EVENT
		  scalar_type empty=G_bar*qvec[g_id+1][rank][t_i][alpha]; 
		  //sum q_sum+=empty;
		  //qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=G_bar*qvec[g_id+1][rank][t_i][alpha];
		  //0.
		  //max
		  if (max_term<empty) 
		    {
		      max_term=empty;
		    }
		  //max		    

		  //sum qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum;
		  if (qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]<max_term)		      
		    {
		      qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]=max_term;
		    }
		  //events within slice rank at time t on alpha virtual branch.
		}
	      if(1)
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
		      //sum scalar_type q_sum=0;
		      //max
		      scalar_type max_term=0;
		      //max

		      //non-leaf directed partition		   
		      if (not is_a_leaf)
			for (int i=0;i<N_parts;i++)
			  {	
			    long int gp_id=gp_ids.at(i);
			    long int gpp_id=gpp_ids.at(i);	    
			    scalar_type pp=p_part.at(i);
			    scalar_type qpe=qvec[gp_id+1][rank][t_i][e];
			    scalar_type qppe=qvec[gpp_id+1][rank][t_i][e];
			    scalar_type Sb_pa_ppe= p_Delta_bar*qvec[gp_id+1][rank][t_i][alpha]*qppe*pp;
			    scalar_type Sb_pe_ppa= p_Delta_bar*qpe*qvec[gpp_id+1][rank][t_i][alpha]*pp;
			    //S_bar EVENT
			    //sum q_sum+= Sb_pa_ppe + Sb_pe_ppa;
			    //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=p_Delta_bar*(qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]+qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]);			  
			    //S_bar.
			    //max
			    if (max_term<Sb_pa_ppe) 
			      {
				max_term=Sb_pa_ppe;
			      }
			    if (max_term<Sb_pe_ppa) 
			      {
				max_term=Sb_pe_ppa;
			      }
			    //max

			    scalar_type D=2*p_delta_e*qpe*qppe*pp;
			    //D EVENT
			    //sum q_sum+= D;
			    //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=2*p_delta_e*qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][e];
			    //D.
			    //max
			    if (max_term<D) 
			      {
				max_term=D;
			      }
			    //max

			  }
		      //sum qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum;
		      if (qvec[g_id+1][tpdt_rank][tpdt_t_i][e]<max_term)		      
			{
			  qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=max_term;
			}

		      scalar_type empty=Get*qvec[g_id+1][rank][t_i][e];
		      //0 EVENT
		      //sum q_sum+=empty;
		      //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=Get*qvec[g_id+1][rank][t_i][e];
		      //0.
		      //max
		      if (max_term<empty) 
			{
			  max_term=empty;
			}
		      //max

		      scalar_type SLb=p_Delta_bar*Eet*qvec[g_id+1][rank][t_i][alpha];
		      //SL_bar EVENT
		      //sum q_sum+=SLb;
		      //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=p_Delta_bar*Eet*qvec[g_id+1][rank][t_i][alpha];
		      //SL_bar.
		      //max
		      if (max_term<SLb) 
			{
			  max_term=SLb;
			}
		      //max

		      //sum qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum;
		      if (qvec[g_id+1][tpdt_rank][tpdt_t_i][e]<max_term)		      
			{
			  qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=max_term;
			}
		      //events within slice rank at time t on branch e. 
		    }	
		 
		}
	      //######################################################################################################################
	      //#########################################INNNER LOOP##################################################################
	      //######################################################################################################################		    
	    }
	}
      gp_ids.clear();
      gpp_ids.clear();
      p_part.clear();
    }
  //del-locs
  g_ids.clear();
  g_id_sizes.clear();
  scalar_type max_term=0;  
  int max_e=-11;
  scalar_type max_t=-11;
  scalar_type max_rank=-11;
  
	
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
	
  
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{
	  //scalar_type t=time_slice_times[rank][t_i];
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {	    
	      int e = time_slices[rank][branch_i];
	      if (max_term<qvec[0][rank][t_i][e]) 
		{
		  max_term=qvec[0][rank][t_i][e];
		  max_e=e;
		  max_t=t_i;
		  max_rank=rank;
		}
	    }
	  if (max_term<qvec[0][rank][t_i][alpha]) 
	    {
	      max_term=qvec[0][rank][t_i][alpha];
	      max_e=alpha;
	      max_t=t_i;
	      max_rank=rank;
	    }	  	  
	}
    } 
  pair<string,scalar_type> return_pair;
  MLRec_events.clear();
  Ttokens.clear();
  register_O(max_e);
  return_pair.first=sample(false,-1,max_t,max_rank,max_e,0,"","",true)+";\n";
  return_pair.second=max_term/root_norm;	
  return return_pair;
}

//used by sample() consider moving to sample.cpp
void exODT_model::register_O(int e)
{
  if (e>-1) branch_counts["count"].at(e)+=1;
  if (e>-1) branch_counts["Os"].at(e)+=1;
}
void exODT_model::register_D(int e)
{
  MLRec_events["D"]+=1;
  if (e>-1) branch_counts["Ds"].at(e)+=1;
}

void exODT_model::register_Tto(int e)
{
  MLRec_events["T"]+=1;
  if (e>-1) branch_counts["Ts"].at(e)+=1; 
}

void exODT_model::register_Tfrom(int e)  
{
  if (e>-1) branch_counts["Tfroms"].at(e)+=1; 
}

void exODT_model::register_L(int e)
{
  MLRec_events["L"]+=1;
  if (e>-1) branch_counts["Ls"].at(e)+=1;
}
void exODT_model::register_S(int e)
{
  MLRec_events["S"]+=1;
  if (e>-1) 
    {
      int f=daughters[e][0];
      int g=daughters[e][1];
      branch_counts["copies"].at(e)+=1;
      branch_counts["count"].at(f)+=1;
      branch_counts["count"].at(g)+=1;  
    }
}
void exODT_model::register_leaf(int e)
{
  if (e>-1) branch_counts["copies"].at(e)+=1;
  //MLRec_events["genes"]+=1;
}

void exODT_model::register_Ttoken(string token)
{
  Ttokens.push_back(token);
}

//ad hoc function should be moved to a future exODT_util.cpp 
void exODT_model::show_counts(string name)
{
  for (map <Node *,int >::iterator it=node_ids.begin();it!=node_ids.end();it++ )
    (*it).first->setBranchProperty("ID",BppString(""));
  
  for (int branch=0;branch<last_branch;branch++)	
    if ( id_nodes.count(branch))
      {
	Node * tmp_node = id_nodes[branch];
	
	stringstream out;
	string old_name = (* (dynamic_cast<const BppString *>(tmp_node->getBranchProperty("ID")))).toSTL();
	//out<< id_ranks[branch];
	if (branch==last_branch-1) out<<"|"<<name<<"|";
	out<<branch_counts[name][branch];
	tmp_node->setBranchProperty("ID",BppString(out.str()));	      
	if (tmp_node->isLeaf())
	  tmp_node->setName(tmp_node->getName()+"_"+out.str());
      }
  cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID") << endl;
  for (map <Node *,int >::iterator it=node_ids.begin();it!=node_ids.end();it++ )
    {
      (*it).first->setBranchProperty("ID",BppString(""));
      if ((*it).first->isLeaf())
	{
	  vector <string> tokens;
	  name=(*it).first->getName();
	  boost::split(tokens,name,boost::is_any_of("_"),boost::token_compress_on);
	  (*it).first->setName(tokens[0]);
	}
    }
  
}

//ad hoc function should be moved to a future exODT_util.cpp 
string exODT_model::counts_string()
{
  stringstream out;
  for (int branch=0;branch<last_branch;branch++)	
    {	
      bool isleaf=false;
      stringstream named_branch;
      if (branch==alpha)
	named_branch<<-1;
      else if (id_ranks[branch]==0)
	{
	  isleaf=true;
	  named_branch<<extant_species[branch];
	}
      else
	named_branch<<id_ranks[branch];
      if (not isleaf)
	out<< "S_internal_branch\t"<< named_branch.str() << "\t" 
	   << branch_counts["Ds"][branch] << "\t"
	   << branch_counts["Ts"][branch] << "\t"
	   << branch_counts["Ls"][branch] << "\t"
	   << branch_counts["copies"][branch] << "\n";
      else
	out<< "S_terminal_branch\t"<< named_branch.str() << "\t" 
	   << branch_counts["Ds"][branch] << "\t"
	   << branch_counts["Ts"][branch] << "\t"
	   << branch_counts["Ls"][branch] << "\t"
	   << branch_counts["copies"][branch] << "\n";
	
    }  
  return out.str();
}

//ad hoc function should be moved to a future exODT_util.cpp 
void exODT_model::show_rates(string name)
{
  for (map <Node *,int >::iterator it=node_ids.begin();it!=node_ids.end();it++ )
    (*it).first->setBranchProperty("ID",BppString(""));
  
  for (int branch=0;branch<last_branch;branch++)	
    if ( id_nodes.count(branch))
      {
	Node * tmp_node = id_nodes[branch];
	
	stringstream out;
	string old_name = (* (dynamic_cast<const BppString *>(tmp_node->getBranchProperty("ID")))).toSTL();
	//out<< id_ranks[branch];
	if (branch==last_branch-1) out<<"|"<<name<<"|";
	if (name=="tau")
	  out<<vector_parameter[name][branch]*vector_parameter["N"][0];
	else
	  out<<vector_parameter[name][branch];
	tmp_node->setBranchProperty("ID",BppString(out.str()));	      
	if (tmp_node->isLeaf())
	  tmp_node->setName(tmp_node->getName()+"_"+out.str());
      }
  
  cout << TreeTemplateTools::treeToParenthesis(*S,false,"ID") << endl;
  for (map <Node *,int >::iterator it=node_ids.begin();it!=node_ids.end();it++ )
    {
      (*it).first->setBranchProperty("ID",BppString(""));
      if ((*it).first->isLeaf())
	{
	  vector <string> tokens;
	  name=(*it).first->getName();
	  boost::split(tokens,name,boost::is_any_of("_"),boost::token_compress_on);
	  (*it).first->setName(tokens[0]);
	}
    }
}
