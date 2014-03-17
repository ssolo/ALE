#include "exODT.h"
using namespace std;
using namespace bpp;

//High memory usage lowmem=false traceback is deprecated!
//The current lowmem=true method uses sample(true) cf. sample.cpp. 
//The general structure of the calculation, and lot of the code, is the same as p(ale) cf. model.cpp.
pair<string,scalar_type> exODT_model::p_MLRec(approx_posterior *ale, bool lowmem)
{
  ale_pointer=ale;
  //cout << "start" << endl;  
  //iterate over directed patitions (i.e. clades) ordered by the number of leaves
  //cout << "start loop" << endl;

  //test  
  //long int tmp_g_id=-1;
  //cout << ale->set2name(ale->id_sets[tmp_g_id]) <<endl;
  //test  
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

  //
  // gene<->species mapping
  for (int i=0;i<(int)g_ids.size();i++)
    {
      long int g_id=g_ids[i];		
      for (int rank=0;rank<last_rank;rank++)
	{
	  int n=time_slices[rank].size();
	  for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	    {
	      scalar_type t=time_slice_times[rank][t_i];
	      for (int branch_i=0;branch_i<n;branch_i++)
		{	    
		  int e = time_slices[rank][branch_i];
		  q[g_id][t][e]=0;
		}
	      q[g_id][t][alpha]=0;
	    }
	}
       
      if (g_id_sizes[i]==1)
	{
        int id = 0;
        boost::dynamic_bitset<> temp = ale->id_sets[g_id];
        for (auto i = 0; i < ale->Gamma_size + 1; ++i) {
           // if ( BipartitionTools::testBit ( temp, i) ) {
            if ( temp[ i ] ) {
                id = i;
                break;
            }
        }
        string gene_name=ale->id_leaves[ id ];
//	  string gene_name=ale->id_leaves[(* (ale->id_sets[g_id].begin()) )];
	  vector <string> tokens;
	  boost::split(tokens,gene_name,boost::is_any_of(string_parameter["gene_name_separators"]),boost::token_compress_on);
	  string species_name;
	  if ((int)scalar_parameter["species_field"]==-1)
	    species_name=tokens[tokens.size()-1];	  
	  else
	    species_name=tokens[(int)scalar_parameter["species_field"]];	  
	  
	    bool found=false;
	    for (map<int,std::string>::iterator tmpit=extant_species.begin();tmpit!=extant_species.end();tmpit++)
	    if ((*tmpit).second==species_name) found = true;	  
	    if (not found) cout << species_name << endl;
	  
	  gid_sps[g_id]=species_name;
	}
    }
  
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
	      scalar_type tpdt,tpdt_nl;
	      //if ( t_i < scalar_parameter["D"]-1 )
	      if ( t_i < (int)time_slice_times[rank].size()-1 )
		tpdt=time_slice_times[rank][t_i+1];
	      else if (rank<last_rank-1)
		tpdt=time_slice_times[rank+1][0];
	      else
		//top of root stem
		tpdt=t_begin[time_slices[rank][0]];

	      if (scalar_parameter["event_node"]==1 and 0)
		tpdt_nl=t;
	      else
		tpdt_nl=tpdt;

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
			    q[g_id][t][e]=1;
			  else
			    q[g_id][t][e]=0;
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
			      //q[g_id][t][e]=0;
			      //max
			      scalar_type max_term=0;
			      step max_step;		       	
			      //max
			    
			      scalar_type SL_fLg=q[g_id][t][f]*Egt;
			      scalar_type SL_Lfg=q[g_id][t][g]*Eft;
			      //SL EVENT
			      // q_sum+=SL_fLg+SL_Lfg;
			      //q[g_id][t][e]=q[g_id][t][f]*Egt + q[g_id][t][g]*Eft;
			      //SL.
			      //max
			      if (max_term<SL_fLg) 
				{
				  max_term= SL_fLg;
				  max_step.e=f;
				  max_step.ep=-1;
				  max_step.epp=-1;
				  max_step.t=t;
				  max_step.rank=rank;
				  max_step.g_id=g_id;
				  max_step.gp_id=-1;
				  max_step.gpp_id=-1;
				  max_step.event="SL";
				}
			      if (max_term<SL_Lfg) 
				{
				  max_term= SL_Lfg;
				  max_step.e=g;
				  max_step.ep=-1;
				  max_step.epp=-1;
				  max_step.t=t;
				  max_step.rank=rank;
				  max_step.g_id=g_id;
				  max_step.gp_id=-1;
				  max_step.gpp_id=-1;
				  max_step.event="SL";
				}
			      //max

			      //non-leaf directed partition
			      if (not is_a_leaf)
				for (int i=0;i<N_parts;i++)
				  {	
				    long int gp_id=gp_ids.at(i);
				    long int gpp_id=gpp_ids.at(i);	    
				    scalar_type pp=p_part.at(i);
				    scalar_type S_pf_ppg=q[gp_id][t][f]*q[gpp_id][t][g]*pp;
				    scalar_type S_ppf_pg=q[gpp_id][t][f]*q[gp_id][t][g]*pp;
				    //S EVENT
				    //q[g_id][t][e]+=q[gp_id][t][f]*q[gpp_id][t][g] +q[gpp_id][t][f]*q[gp_id][t][g];
				    //sum q_sum+= S_pf_ppg + S_ppf_pg;
				    //S.
				    //max
				    if (max_term<S_pf_ppg) 
				      {
					max_term= S_pf_ppg;
					max_step.e=-1;
					max_step.ep=f;
					max_step.epp=g;
					max_step.t=t;
					max_step.rank=rank;
					max_step.g_id=-1;
					max_step.gp_id=gp_id;
					max_step.gpp_id=gpp_id;
					max_step.event="S";
				      }
				    if (max_term<S_ppf_pg) 
				      {
					max_term=S_ppf_pg;
					max_step.e=-1;
					max_step.ep=g;
					max_step.epp=f;
					max_step.t=t;
					max_step.rank=rank;
					max_step.g_id=-1;
					max_step.gp_id=gp_id;
					max_step.gpp_id=gpp_id;
					max_step.event="S";
				      }
				    //max

				  }

			      //sum q[g_id][t][e]=q_sum; 
			      q[g_id][t][e]=max_term; 
			      if (not lowmem) q_step[g_id][t][e]=max_step; 
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
	      if(1)
		{
		  //boundaries for branch alpha virtual branch  
		  //boundary at present
		  if (t==0)
		    q[g_id][t][alpha]=0;
		  //boundary between slice rank and rank-1 slice is trivial	
		  ;//q[g_id][t][alpha]=q[g_id][t][alpha];	  
		  //boundaries for branch alpha virtual branch.  

		  //events within slice rank at time t on alpha virtual branch
		  scalar_type G_bar=Ge[-1][t];
		  q[g_id][tpdt][alpha]=0;
		  //sum scalar_type q_sum=0;
		  //sum scalar_type q_sum_nl=0;
		  //max
		  scalar_type max_term=0;
		  scalar_type max_term_nl=0;
		  step max_step;
		  step max_step_nl;
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
			    scalar_type T_ep_app=p_Ntau_e*q[gp_id][t][e]*q[gpp_id][t][alpha]*pp;
			    scalar_type T_ap_epp=p_Ntau_e*q[gp_id][t][alpha]*q[gpp_id][t][e]*pp;
			    //Tb EVENT
			    //sum q_sum_nl+=T_ep_app+T_ap_epp;
			    //q[g_id][tpdt][alpha]+=p_Ntau_e*(q[gp_id][t][e]*q[gpp_id][t][alpha]+q[gp_id][t][alpha]*q[gpp_id][t][e]);
			    //Tb.
			    //max
			    if (max_term_nl<T_ep_app) 
			      {
				max_term_nl=T_ep_app;
				max_step_nl.e=-1;
				max_step_nl.ep=e;
				max_step_nl.epp=alpha;
				max_step_nl.t=t;
				max_step_nl.rank=rank;
				max_step_nl.g_id=-1;
				max_step_nl.gp_id=gp_id;
				max_step_nl.gpp_id=gpp_id;
				max_step_nl.event="Tb";			      
			      }
			    if (max_term_nl<T_ap_epp) 
			      {
				max_term_nl=T_ap_epp;
				max_step_nl.e=-1;
				max_step_nl.ep=alpha;
				max_step_nl.epp=e;
				max_step_nl.t=t;
				max_step_nl.rank=rank;
				max_step_nl.g_id=-1;
				max_step_nl.gp_id=gp_id;
				max_step_nl.gpp_id=gpp_id;
				max_step_nl.event="Tb";			      
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
			scalar_type Sb=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha])*pp;
			//S_bar EVENT
			//sum q_sum_nl+=Sb;
			//q[g_id][tpdt][alpha]+=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha]);
			//S_bar.
			//max
			if (max_term_nl<Sb) 
			  {
			    max_term_nl=Sb;
			    max_step_nl.e=-1;
			    max_step_nl.ep=alpha;
			    max_step_nl.epp=alpha;
			    max_step_nl.t=t;
			    max_step_nl.rank=rank;
			    max_step_nl.g_id=-1;
			    max_step_nl.gp_id=gp_id;
			    max_step_nl.gpp_id=gpp_id;
			    max_step_nl.event="Sb";			      
			  }
			//max

		      }	    
	      
		  for (int branch_i=0;branch_i<n;branch_i++)			  
		    {
		      int e = time_slices[rank][branch_i];		
		      scalar_type tau_e=vector_parameter["tau"][e];
		      scalar_type p_Ntau_e=tau_e*Delta_t;
		      scalar_type TLb=p_Ntau_e*Ebar*q[g_id][t][e];
		      //TL_bar EVENT
		      //sum q_sum+=TLb;
		      //q[g_id][tpdt][alpha]+=p_Ntau_e*Ebar*q[g_id][t][e];
		      //TL_bar.
		      //max
		      if (max_term_nl<TLb) 
			{
			  max_term_nl=TLb;
			  max_step_nl.e=e;
			  max_step_nl.ep=-1;
			  max_step_nl.epp=-1;
			  max_step_nl.t=t;
			  max_step_nl.rank=rank;
			  max_step_nl.g_id=g_id;
			  max_step_nl.gp_id=-1;
			  max_step_nl.gpp_id=-1;
			  max_step_nl.event="TLb";			      			
			}
		      //max		    
		    }
	      
		  //sum q[g_id][tpdt_nl][alpha]+=q_sum_nl;
		  if (q[g_id][tpdt_nl][alpha]<max_term_nl)		      
		    {
		      q[g_id][tpdt_nl][alpha]=max_term_nl;
		      if (not lowmem) q_step[g_id][tpdt_nl][alpha]=max_step_nl;
		    }
	      
		  //0 EVENT
		  scalar_type empty=G_bar*q[g_id][t][alpha]; 
		  //sum q_sum+=empty;
		  //q[g_id][tpdt][alpha]+=G_bar*q[g_id][t][alpha];
		  //0.
		  //max
		  if (max_term<empty) 
		    {
		      max_term=empty;
		      max_step.e=alpha;
		      max_step.ep=-1;
		      max_step.epp=-1;
		      max_step.t=t;
		      max_step.rank=rank;
		      max_step.g_id=g_id;
		      max_step.gp_id=-1;
		      max_step.gpp_id=-1;
		      max_step.event="0";			      
		    }
		  //max		    

		  //sum q[g_id][tpdt][alpha]+=q_sum;
		  if (q[g_id][tpdt][alpha]<max_term)		      
		    {
		      q[g_id][tpdt][alpha]=max_term;
		      if (not lowmem) q_step[g_id][tpdt][alpha]=max_step;
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
		      q[g_id][tpdt][e]=0;
		      //sum scalar_type q_sum=0;
		      //sum scalar_type q_sum_nl=0;
		      //max
		      scalar_type max_term=0;
		      scalar_type max_term_nl=0;
		      step max_step;
		      step max_step_nl;
		      //max

		      //non-leaf directed partition		   
		      if (not is_a_leaf)
			for (int i=0;i<N_parts;i++)
			  {	
			    long int gp_id=gp_ids.at(i);
			    long int gpp_id=gpp_ids.at(i);	    
			    scalar_type pp=p_part.at(i);
			    scalar_type qpe=q[gp_id][t][e];
			    scalar_type qppe=q[gpp_id][t][e];
			    scalar_type Sb_pa_ppe= p_Delta_bar*q[gp_id][t][alpha]*qppe*pp;
			    scalar_type Sb_pe_ppa= p_Delta_bar*qpe*q[gpp_id][t][alpha]*pp;
			    //S_bar EVENT
			    //sum q_sum_nl+= Sb_pa_ppe + Sb_pe_ppa;
			    //q[g_id][tpdt][e]+=p_Delta_bar*(q[gp_id][t][alpha]*q[gpp_id][t][e]+q[gp_id][t][e]*q[gpp_id][t][alpha]);			  
			    //S_bar.
			    //max
			    if (max_term_nl<Sb_pa_ppe) 
			      {
				max_term_nl=Sb_pa_ppe;
				max_step_nl.e=-1;
				max_step_nl.ep=alpha;
				max_step_nl.epp=e;
				max_step_nl.t=t;
				max_step_nl.rank=rank;
				max_step_nl.g_id=-1;
				max_step_nl.gp_id=gp_id;
				max_step_nl.gpp_id=gpp_id;
				max_step_nl.event="Sb";
			      }
			    if (max_term_nl<Sb_pe_ppa) 
			      {
				max_term_nl=Sb_pe_ppa;
				max_step_nl.e=-1;
				max_step_nl.ep=e;
				max_step_nl.epp=alpha;
				max_step_nl.t=t;
				max_step_nl.rank=rank;
				max_step_nl.g_id=-1;
				max_step_nl.gp_id=gp_id;
				max_step_nl.gpp_id=gpp_id;
				max_step_nl.event="Sb";
			      }
			    //max

			    scalar_type D=2*p_delta_e*qpe*qppe*pp;
			    //D EVENT
			    //sum q_sum_nl+= D;
			    //q[g_id][tpdt][e]+=p_delta_e*q[gp_id][t][e]*q[gpp_id][t][e];
			    //D.
			    //max
			    if (max_term_nl<D) 
			      {
				max_term_nl=D;
				max_step_nl.e=-1;
				max_step_nl.ep=e;
				max_step_nl.epp=e;
				max_step_nl.t=t;
				max_step_nl.rank=rank;
				max_step_nl.g_id=-1;
				max_step_nl.gp_id=gp_id;
				max_step_nl.gpp_id=gpp_id;
				max_step_nl.event="D";
			      }
			    //max

			  }
		      //sum q[g_id][tpdt_nl][e]+=q_sum_nl;
		      if (q[g_id][tpdt_nl][e]<max_term_nl)		      
			{
			  q[g_id][tpdt_nl][e]=max_term_nl;
			  if (not lowmem) q_step[g_id][tpdt_nl][e]=max_step_nl;
			}

		      scalar_type empty=Get*q[g_id][t][e];
		      //0 EVENT
		      //sum q_sum+=empty;
		      //q[g_id][tpdt][e]=Get*q[g_id][t][e];
		      //0.
		      //max
		      if (max_term<empty) 
			{
			  max_term=empty;
			  max_step.e=e;
			  max_step.ep=-1;
			  max_step.epp=-1;
			  max_step.t=t;
			  max_step.rank=rank;
			  max_step.g_id=g_id;
			  max_step.gp_id=-1;
			  max_step.gpp_id=-1;
			  max_step.event="0";			      
			}
		      //max

		      scalar_type SLb=p_Delta_bar*Eet*q[g_id][t][alpha];
		      //SL_bar EVENT
		      //sum q_sum+=SLb;
		      //q[g_id][tpdt][e]+=p_Delta_bar*Eet*q[g_id][t][alpha];
		      //SL_bar.
		      //max
		      if (max_term<SLb) 
			{
			  max_term=SLb;
			  max_step.e=alpha;
			  max_step.ep=-1;
			  max_step.epp=-1;
			  max_step.t=t;
			  max_step.rank=rank;
			  max_step.g_id=g_id;
			  max_step.gp_id=-1;
			  max_step.gpp_id=-1;
			  max_step.event="SLb";			      
			}
		      //max

		      //sum q[g_id][tpdt][e]+=q_sum;
		      if (q[g_id][tpdt][e]<max_term)		      
			{
			  q[g_id][tpdt][e]=max_term;
			  if (not lowmem) q_step[g_id][tpdt][e]=max_step;
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
  if (not lowmem)   
    {
      return traceback();	
    }
  else 
    {
      //cout << "LOWMEM" <<endl;
      scalar_type max_term=0;  
      long int g_id=-1;
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
	      scalar_type t=time_slice_times[rank][t_i];
	      for (int branch_i=0;branch_i<n;branch_i++)
		{	    
		  int e = time_slices[rank][branch_i];
		  if (max_term<q[g_id][t][e]) 
		    {
		      max_term=q[g_id][t][e];
		      max_e=e;
		      max_t=t_i;
		      max_rank=rank;
		    }
		}
	      if (max_term<q[g_id][t][alpha]) 
		{
		  max_term=q[g_id][t][alpha];
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
}

//deprecated
pair<string,scalar_type> exODT_model::traceback()
{
  stringstream signal_stream;
  scalar_type max_term=0;  
  long int g_id=-1;
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
	  scalar_type t=time_slice_times[rank][t_i];
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {	    
	      int e = time_slices[rank][branch_i];
	      signal_stream << g_id<<" "<<t<<" "<<e<<" "<<q[g_id][t][e]<<endl;
	      if (max_term<q[g_id][t][e]) 
		{
		  max_term=q[g_id][t][e];
		  max_e=e;
		  max_t=t;
		  max_rank=rank;
		}
	    }
	  signal_stream << g_id<<" "<<t<<" "<<alpha<<" "<<q[g_id][t][alpha]<<endl;
	  if (max_term<q[g_id][t][alpha]) 
	    {
	      max_term=q[g_id][t][alpha];
	      max_e=alpha;
	      max_t=t;
	      max_rank=rank;
	    }	  	  
	}
    } 
  signal_string=signal_stream.str();
  pair<string,scalar_type> return_pair;
  MLRec_events.clear();
  Ttokens.clear();
  register_O(max_e);
  return_pair.first=traceback(g_id,max_t,max_rank,max_e,0,"")+";\n";
  return_pair.second=max_term/root_norm;

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
  return return_pair;
}

//deprecated
string exODT_model::traceback(long int g_id,scalar_type t,scalar_type rank,int e,scalar_type branch_length,string branch_events, string transfer_token)
{
  /*
  if (e==alpha)
    cout << "b "<<-1;
  else if (id_ranks[e]==0)
    cout << "b "<<extant_species[e];
  else
    cout << "b "<< id_ranks[e];
  
  cout << "from " << e << " " << t << " " << endl; 
  cout << "g_id "<< g_id <<" is " <<ale_pointer->set2name(ale_pointer->id_sets[g_id]) << endl;
  */
  step max_step=q_step[g_id][t][e];

  //cout << max_step.event<<endl;
  
  //cout << t << "==0 and " << (int)( tmpset.size() ) << "==1 and " <<e <<" !=-1 "<<endl;
  scalar_type new_branch_length=branch_length+t-max_step.t;
  if (max_step.t!=max_step.t)
    new_branch_length=branch_length;
    int size = 0;
    boost::dynamic_bitset<>  temp = ale_pointer->id_sets[g_id];
    for (auto i = 0; i < ale_pointer->Gamma_size + 1; ++i) {
      //  if ( BipartitionTools::testBit ( temp, i) ) {
        if ( temp[ i] ) {
        size++;
        }
    }
  if (t==0 and size==1 and e!=-1)
    {
      register_leaf(e);
      stringstream branch_string;
      if (scalar_parameter["leaf_events"]==1) branch_string<<branch_events;
      branch_string <<":"<<new_branch_length;
      return ale_pointer->set2name(ale_pointer->id_sets[g_id])+branch_string.str();
    }
  /*
  if (ale_pointer->Bip_counts[g_id]>0)
    new_branch_length=ale_pointer->Bip_bls[g_id]/ale_pointer->Bip_counts[g_id];
  else
    new_branch_length=ale_pointer->Bip_bls[g_id]/ale_pointer->observations;
  */
  if (max_step.event=="D" or max_step.event=="Tb" or max_step.event=="S" or max_step.event=="Sb")
    {

      //cout << max_step.event << " " << max_step.ep << "-" << max_step.epp <<" "<<max_step.t <<endl; 
      stringstream transfer_token_stream;
      transfer_token_stream<<"";
      stringstream branch_string;
      if (max_step.event=="S")
	{
	  register_S(e);
	  branch_string<< branch_events//<<max_step.event<<"@"
		       <<"."<<id_ranks[e]<<":"<<max(new_branch_length,(scalar_type)0.0); 
	}
      else
	{
	  if (max_step.event=="Tb")
	    {
	      //cout << "T";
	      int this_e;
	      if (max_step.ep==alpha)
		this_e=max_step.epp;
	      else
		this_e=max_step.ep;
	      stringstream named_branch;	      
	      if (this_e==alpha)
		named_branch<<-1;
	      else if (id_ranks[this_e]==0)
		named_branch<<extant_species[this_e];
	      else
		named_branch<<id_ranks[this_e];
	      // Tto
	      register_Tto(this_e);
	      stringstream tmp;
	      tmp<<max_step.rank<<"|"<<t<<"|"<<named_branch.str();
	      //tmp<<max_step.t<<"|"<<this_e;
	      register_Ttoken(transfer_token+"|"+tmp.str());
	      // Tto
	      
	      branch_string<< branch_events<<max_step.event<<"@"<<max_step.rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0); 
	    }
	  else if (max_step.event=="Sb")
	    {
	      int this_e;
	      if (max_step.ep==alpha)
		this_e=max_step.epp;
	      else
		this_e=max_step.ep;
	      // Tfrom
	      register_Tfrom(this_e);
	      // Tfrom
	      stringstream named_branch;
	      if (this_e==alpha)
		named_branch<<-1;
	      else if (id_ranks[this_e]==0)
		named_branch<<extant_species[this_e];
	      else
		named_branch<<id_ranks[this_e];
	      
	      transfer_token_stream<< rank<<"|"<<t<<"|"<<named_branch.str();
	      //transfer_token_stream<< t<<"|"<<this_e;

	      branch_string<< branch_events<<"T@"<<rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0); 	    
	    }
	  else
	    {
	      //cout << "D";
	      register_D(e);	      
	      stringstream named_branch;
	      if (e==alpha)
		named_branch<<-1;
	      else if (id_ranks[e]==0)
		named_branch<<extant_species[e];
	      else
		named_branch<<id_ranks[e];
	    
	      branch_string<< branch_events<<max_step.event<<"@"<<rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0); 
	    }
	}
      if ( transfer_token_stream.str()=="")
	return "("+
	  traceback(max_step.gp_id,max_step.t,max_step.rank,max_step.ep,0,"",transfer_token)
	  +","+
	  traceback(max_step.gpp_id,max_step.t,max_step.rank,max_step.epp,0,"",transfer_token)
	  +")"+branch_string.str();
      else
	if(max_step.ep==alpha)
	  return "("+
	    traceback(max_step.gp_id,max_step.t,max_step.rank,max_step.ep,0,"",transfer_token_stream.str())
	    +","+
	    traceback(max_step.gpp_id,max_step.t,max_step.rank,max_step.epp,0,"",transfer_token)
	    +")"+branch_string.str();
	else
	  return "("+
	    traceback(max_step.gp_id,max_step.t,max_step.rank,max_step.ep,0,"",transfer_token)
	    +","+
	    traceback(max_step.gpp_id,max_step.t,max_step.rank,max_step.epp,0,"",transfer_token_stream.str())
	    +")"+branch_string.str();
	  
    }
  else if ( max_step.event=="TLb" or max_step.event=="SL" or max_step.event=="SLb" or max_step.event=="0")
    {
      //cout << max_step.event << " " << max_step.e <<" "<<max_step.t <<endl; 

      stringstream branch_string;
      stringstream transfer_token_stream;
      transfer_token_stream <<"";
      
      branch_string<< branch_events;
      if (max_step.event!="0")
	{
	  if (max_step.event=="SL")
	    {
	      register_S(e);	      
	      int f=daughters[e][0];
	      int g=daughters[e][1];
	      if (max_step.e==f)
		register_L(g);
	      else
		register_L(f);
	      branch_string<<"."//<<max_step.event<<"@"
			   <<id_ranks[e];
	    }
	  else
	    {
	      if (max_step.event=="TLb")
		{
		  register_Tto(max_step.e);
		  stringstream tmp;

		  stringstream named_branch;
		  if (max_step.e==alpha)
		    named_branch<<-1;
		  else if (id_ranks[max_step.e]==0)
		    named_branch<<extant_species[max_step.e];
		  else
		    named_branch<<id_ranks[max_step.e];

		  tmp<<max_step.rank<<"|"<< t << "|"<< named_branch.str();
		  //tmp<<max_step.t<<"|"<<max_step.e;
		  register_Ttoken(transfer_token+"|"+tmp.str());
		  transfer_token="";
	      
		  branch_string<<""//<<max_step.event
			       <<"@"<<max_step.rank<<"|"<<named_branch.str();
		}
	      else  if (max_step.event=="SLb")
		{
		  //cout << "L";
		  register_L(e);
		  register_Tfrom(e);
		  stringstream named_branch;
		  if (e==alpha)
		    named_branch<<-1;
		  else if (id_ranks[e]==0)
		    named_branch<<extant_species[e];
		  else
		    named_branch<<id_ranks[e];

		  transfer_token_stream<< rank<<"|"<<t<<"|"<<named_branch.str();
		  //transfer_token_stream<< max_step.t<<"|"<<e;

		  branch_string<<".T"//<<max_step.event
			       <<"@"<<rank<<"|"<<named_branch.str();
		}
	    }
	}
      if (transfer_token_stream.str()=="")
	return traceback(max_step.g_id, max_step.t,max_step.rank, max_step.e, new_branch_length, branch_string.str(),transfer_token);
      else
	return traceback(max_step.g_id, max_step.t,max_step.rank, max_step.e, new_branch_length, branch_string.str(),transfer_token_stream.str());
    }
  else
    {
      cout <<  "me " <<max_step.event << endl;
      cout << "error "  <<endl;
      cout << " g_id " << g_id;
      cout << " t " << t;
      cout << " e " << e; 
      cout << " l " << branch_length;
      cout << " str " << branch_events;
      cout << endl;
      cout << ale_pointer->constructor_string <<endl;
      signal=-11;
    }
  return "error";
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
