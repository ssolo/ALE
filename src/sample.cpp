#include "exODT.h"
using namespace std;
using namespace bpp;
//
//consider reimplemtation for clarity!
//
//The general structure of the calculation, and lot of the code, is the same as p(ale) cf. model.cpp.
//(this could be made more clear)
string exODT_model::sample(bool max_rec)
{
  MLRec_events.clear();
  Ttokens.clear();
  
  //scalar_type beta=1;
  scalar_type root_resum=0;
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{
	  scalar_type t=time_slice_times[rank][t_i];		
	  
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {
	      int e = time_slices[rank][branch_i];
	      root_resum+=q[-1][t][e];
	    }	     
	  root_resum+=q[-1][t][alpha]; 
	}
    }
  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  scalar_type root_reresum=0;
  int root_e=-1;
  int root_rank=-1;
  int root_t_i=-1;

  int max_root_e=-1;
  int max_root_rank=-1;
  int max_root_t_i=-1;

  scalar_type max_resum=0;
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{
	  scalar_type t=time_slice_times[rank][t_i];		
	  
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {
	      int e = time_slices[rank][branch_i];
	      root_reresum+=q[-1][t][e];	      
	      if (max_resum<q[-1][t][e])
		{
		  max_resum=q[-1][t][e];
		  max_root_e=e;
		  max_root_rank=rank;
		  max_root_t_i=t_i;		  
		}
	      if (r*root_resum<root_reresum and root_t_i==-1)
		{
		  root_e=e;
		  root_rank=rank;
		  root_t_i=t_i;
		}
	    }	     
	  root_reresum+=q[-1][t][alpha];
	  if (max_resum<q[-1][t][alpha])
	    {
	      max_resum=q[-1][t][alpha];
	      max_root_e=alpha;
	      max_root_rank=rank;
	      max_root_t_i=t_i;		  
	    }
	  if (r*root_resum<root_reresum and root_t_i==-1)
	    {
	      root_e=alpha;
	      root_rank=rank;
	      root_t_i=t_i;
	    }	  
	}
    }
  if (max_rec)
    {
      root_e=max_root_e;
      root_rank=max_root_rank;
      root_t_i=max_root_t_i;      
    }
  register_O(root_e);
  return sample(false,-1,root_t_i,root_rank,root_e,0,"","",max_rec)+";";
  //del-locs
}

//
//consider reimplemtation for clarity!
//
//The general structure of the calculation, and lot of the code, is the same as p(ale) cf. model.cpp.
//(this could be made more clear)
string exODT_model::sample(bool S_node,long int g_id,int t_i,scalar_type rank,int e,scalar_type branch_length,string branch_events, string transfer_token,bool max_rec)
{
  // it could be nice to implemant a sampling temperature ?  
  //scalar_type beta=1;
  stringstream topptmp;
  if (e==alpha)
    topptmp<<-1;
  else if (id_ranks[e]==0)
    topptmp<<extant_species[e];
  else
    topptmp<<id_ranks[e];

  approx_posterior * ale=ale_pointer;
  bool is_a_leaf=false; 

  int size = 0;
  boost::dynamic_bitset<> temp;
  //  if (g_id!=-1) { //We are not at the root bipartition
  temp = ale->id_sets.at( g_id );
  for (int i = 0; i < ale->Gamma_size + 1; ++i) {
    // if ( BipartitionTools::testBit ( temp, i) ) {
    if ( temp[ i ] ) {
      size++;
    }
  }
        
        
  // if ((int)(ale->id_sets[g_id].size())==1)
  if (size == 1)
    is_a_leaf=true;
  //  }
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
      //root bipartition needs to be handled seperatly
      map<set<long int>,int> bip_parts;
      for (map <long int,scalar_type> :: iterator it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
	{
          long int gp_id=(*it).first;
          boost::dynamic_bitset<> gamma = ale->id_sets[gp_id];
          
          boost::dynamic_bitset<> not_gamma = ~gamma;
          not_gamma[0] = 0;

          /* for (auto i = 0; i < ale->nbint; ++i) {
	     not_gamma[i] = 0;
	     }
	     BipartitionTools::bitNot(not_gamma, gamma, ale->nbint);*/
          /*
	    for (set<int>::iterator st=ale->Gamma.begin();st!=ale->Gamma.end();st++)
	    if (gamma.count(*st)==0)
	    not_gamma.insert(*st);*/
          long int gpp_id = ale->set_ids[not_gamma];
          set <long int> parts;
          parts.insert(gp_id);
          parts.insert(gpp_id);
          bip_parts[parts]=1;
          //  gamma.clear();
          //  not_gamma.clear();
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
  int n=time_slices[rank].size();
  //######################################################################################################################
  //######################################### INNER LOOP #################################################################
  //######################################################################################################################
      
  vector <step> sample_steps;
  vector <scalar_type> sample_ps;
  scalar_type resum=0;
  scalar_type t;

  //scalar_type t_to=time_slice_times[rank][t_i];

  //int rank_to=rank;
  //int t_i_to=t_i; 
  bool set_S_node=false;
  // proceed a single "D" subslice
  if(t_i>0)
    {
      rank=rank;
      t_i-=1;
    }
  // at boundaries  
  else if (rank>0)
    {
      if (S_node)
	{
	  ;    
	}
      //if e defines the time slice we have to look at speciaitons 
      else if (e==time_slices[rank][n-1])
	{
	  set_S_node=true;
	}
      else
	{
	  rank-=1;
	  t_i=time_slice_times[rank].size()-1;
	}
    }
  else
    {
      rank=-1;      
      t_i=-1;      
      if (is_a_leaf && extant_species[e]==gid_sps[g_id])			  
	{
	  resum=1;
	  sample_ps.push_back(1);
	  step step;
	  step.e=e;
	  step.ep=-1;
	  step.epp=-1;
	  step.t=0;
	  step.rank=0;
	  step.g_id=g_id;
	  step.gp_id=-1;
	  step.gpp_id=-1;
	  step.event="0";			
	  sample_steps.push_back(step);      
	}//q[g_id][t][e]=1;

    }
  if (rank>-1)
    {

      t=time_slice_times[rank][t_i];
      scalar_type tpdt;//,tpdt_nl;      
      if ( t_i < (int)time_slice_times[rank].size()-1 )
	tpdt=time_slice_times[rank][t_i+1];
      else if (rank<last_rank-1)
	tpdt=time_slice_times[rank+1][0];
      else
	//top of root stem
	tpdt=t_begin[time_slices[rank][0]];

      /*
	if (scalar_parameter["event_node"]==1 and 0)
	tpdt_nl=t;
	else
	tpdt_nl=tpdt;
      */
      //root
      scalar_type Delta_t=tpdt-t;
      //scalar_type N=vector_parameter["N"][rank];
      scalar_type Delta_bar=vector_parameter["Delta_bar"][rank];
      //scalar_type Lambda_bar=vector_parameter["Lambda_bar"][rank];
      scalar_type p_Delta_bar=Delta_bar*Delta_t;			     
      scalar_type Ebar=Ee[-1][t];

      if(e==alpha)
	{
	  //boundaries for branch alpha virtual branch  
	  //boundary at present
	  if (t==0)
	    {
	      resum+=0;	      
	      sample_ps.push_back(0);
	      step step;
	      step.e=alpha;
	      step.ep=-1;
	      step.epp=-1;
	      step.t=t;
	      step.rank=rank;
	      step.g_id=g_id;
	      step.gp_id=-1;
	      step.gpp_id=-1;
	      step.event="0";			      
	      sample_steps.push_back(step);
	    }//q[g_id][t][alpha]=0;
	  //boundary between slice rank and rank-1 slice is trivial
	  //trivial
	  if (S_node )//and 0!?
	    {
	      resum=1;
	      if(1)
		{
		  sample_ps.push_back(1);
		  step step;
		  step.e=alpha;
		  step.ep=-1;
		  step.epp=-1;
		  step.t=t;
		  step.rank=rank;
		  step.g_id=g_id;
		  step.gp_id=-1;
		  step.gpp_id=-1;
		  step.event="0";			      
		  sample_steps.push_back(step);
		}
	      ;//q[g_id][t][e]=q[g_id][t][e];
	    }
	  //q[g_id][t][alpha]=q[g_id][t][alpha];	  
	  else
	    {
	      //cout << " here " <<endl;
	      //boundaries for branch alpha virtual branch.  
	      //events within slice rank at time t on alpha virtual branch
	      scalar_type G_bar=Ge[-1][t];
	      //q[g_id][tpdt][alpha]=0;
	      //scalar_type q_sum=0;
	      //scalar_type q_sum_nl=0;
	      for (int branch_i=0;branch_i<n;branch_i++)			  
		{
		  int e = time_slices[rank][branch_i];		
		  scalar_type tau_e=vector_parameter["tau"][e];		    
		  scalar_type p_Ntau_e=tau_e*Delta_t;

		  //non-leaf directed partition
		  if (not is_a_leaf)
		    for (int i=0;i<N_parts;i++)
		      {	
			long int gp_id=gp_ids[i];
			long int gpp_id=gpp_ids[i];	    
			scalar_type pp=p_part[i];				    
			scalar_type T_ep_app=p_Ntau_e*q[gp_id][t][e]*q[gpp_id][t][alpha]*pp;
			scalar_type T_ap_epp=p_Ntau_e*q[gp_id][t][alpha]*q[gpp_id][t][e]*pp;
			//T EVENT
			resum+=T_ep_app;
			if (1)
			  {
			    sample_ps.push_back(T_ep_app);
			    step step;
			    step.e=-1;
			    step.ep=e;
			    step.epp=alpha;
			    step.t=t;
			    step.rank=rank;
			    step.g_id=-1;
			    step.gp_id=gp_id;
			    step.gpp_id=gpp_id;
			    step.event="Tb";
			    sample_steps.push_back(step);	       	      
			  }

			resum+=T_ap_epp;
			if (1)
			  {
			    sample_ps.push_back(T_ap_epp);
			    step step;
			    step.e=-1;
			    step.ep=alpha;
			    step.epp=e;
			    step.t=t;
			    step.rank=rank;
			    step.g_id=-1;
			    step.gp_id=gp_id;
			    step.gpp_id=gpp_id;
			    step.event="Tb";			      
			    sample_steps.push_back(step);
			  }

			//q_sum_nl+=
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
		    scalar_type Sb=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha])*pp;
		    //S_bar EVENT
		    resum+=Sb;
		    if (1)
		      {
			sample_ps.push_back(Sb);
			step step;
			step.e=-1;
			step.ep=alpha;
			step.epp=alpha;
			step.t=t;
			step.rank=rank;
			step.g_id=-1;
			step.gp_id=gp_id;
			step.gpp_id=gpp_id;
			step.event="Sb";			      
			sample_steps.push_back(step);
		      }

		    //q_sum_nl+=Sb;
		    //q[g_id][tpdt][alpha]+=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha]);
		    //S_bar.

		  }	    
		
	      //q[g_id][tpdt_nl][alpha]+=q_sum_nl;
	  
	      for (int branch_i=0;branch_i<n;branch_i++)			  
		{
		  int e = time_slices[rank][branch_i];		
		  scalar_type tau_e=vector_parameter["tau"][e];
		  scalar_type p_Ntau_e=tau_e*Delta_t;
		  scalar_type TLb=p_Ntau_e*Ebar*q[g_id][t][e];
		  //TL_bar EVENT
		  resum+=TLb;
		  if (1)
		    {
		      sample_ps.push_back(TLb);
		      step step;
		      step.e=e;
		      step.ep=-1;
		      step.epp=-1;
		      step.t=t;
		      step.rank=rank;
		      step.g_id=g_id;
		      step.gp_id=-1;
		      step.gpp_id=-1;
		      step.event="TLb";			      
		      sample_steps.push_back(step);
		    }

		  //q_sum+=TLb;
		  //q[g_id][tpdt][alpha]+=p_Ntau_e*Ebar*q[g_id][t][e];
		  //TL_bar.

		}

	      //0 EVENT
	      scalar_type empty=G_bar*q[g_id][t][alpha]; 
	      resum+=empty;
	      if (1)
		{
		  sample_ps.push_back(empty);
		  step step;
		  step.e=alpha;
		  step.ep=-1;
		  step.epp=-1;
		  step.t=t;
		  step.rank=rank;
		  step.g_id=g_id;
		  step.gp_id=-1;
		  step.gpp_id=-1;
		  step.event="0";		
		  sample_steps.push_back(step);	      
		}

	      //q_sum+=empty;

	      //q[g_id][tpdt][alpha]+=G_bar*q[g_id][t][alpha];
	      //0.

	      //q[g_id][tpdt][alpha]+=q_sum;
	      //events within slice rank at time t on alpha virtual branch.
	    }
	}
      else
	{

	  //int e = time_slices[rank][branch_i];
	  scalar_type Get=Ge[e][t];
	  scalar_type Eet=Ee[e][t];	
	  scalar_type delta_e=vector_parameter["delta"][e];
	  scalar_type p_delta_e=delta_e*Delta_t;
	  if (S_node)
	    {
	      //boundaries for branch e
	      //boundary at present
	      if (t==0)
		{
		  ;
		}
	      //boundary between slice rank and rank-1
	      else if (t_i==0 )
		{
		  //terminating branch is last in time_slices and defines a represented speciation 
		  if (e==time_slices[rank][n-1] && rank>0)
		    {
		      int f=daughters[e][0];
		      int g=daughters[e][1];
		      scalar_type Eft=Ee[f][t];
		      scalar_type Egt=Ee[g][t];
		      
		      //scalar_type q_sum=0;
		      //q[g_id][t][e]=0;
		      
		      scalar_type SL_fLg=q[g_id][t][f]*Egt;
		      scalar_type SL_Lfg=q[g_id][t][g]*Eft;
		      //SL EVENT
		      resum+=SL_fLg;
		      if(1)
			{
			  sample_ps.push_back(SL_fLg);
			  step step;
			  step.e=f;
			  step.ep=-1;
			  step.epp=-1;
			  step.t=t;
			  step.rank=rank;
			  step.g_id=g_id;
			  step.gp_id=-1;
			  step.gpp_id=-1;
			  step.event="SL";			      
			  sample_steps.push_back(step);
			}

		      resum+=SL_Lfg;
		      if(1)
			{
			  sample_ps.push_back(SL_Lfg);
			  step step;
			  step.e=g;
			  step.ep=-1;
			  step.epp=-1;
			  step.t=t;
			  step.rank=rank;
			  step.g_id=g_id;
			  step.gp_id=-1;
			  step.gpp_id=-1;
			  step.event="SL";		
			  sample_steps.push_back(step);	      
			}

		      //q_sum+=SL_fLg+SL_Lfg;
		      //q[g_id][t][e]=q[g_id][t][f]*Egt + q[g_id][t][g]*Eft;
		      //SL.

		      //non-leaf directed partition
		      if (not is_a_leaf)
			for (int i=0;i<N_parts;i++)
			  {	
			   
			    long int gp_id=gp_ids[i];
			    long int gpp_id=gpp_ids[i];	    
			    scalar_type pp=p_part[i];
			    scalar_type S_pf_ppg=q[gp_id][t][f]*q[gpp_id][t][g]*pp;
			    scalar_type S_ppf_pg=q[gpp_id][t][f]*q[gp_id][t][g]*pp;
			    //S EVENT
			    //q[g_id][t][e]+=q[gp_id][t][f]*q[gpp_id][t][g] +q[gpp_id][t][f]*q[gp_id][t][g];
			    resum+=S_pf_ppg;
			    if(1)
			      {
				sample_ps.push_back(S_pf_ppg);
				step step;
				step.e=-1;
				step.ep=f;
				step.epp=g;
				step.t=t;
				step.rank=rank;
				step.g_id=-1;
				step.gp_id=gp_id;
				step.gpp_id=gpp_id;
				step.event="S";			      
				sample_steps.push_back(step);
			      }


			    resum+=S_ppf_pg;
			    if(1)
			      {
				sample_ps.push_back(S_ppf_pg);
				step step;
				step.e=-1;
				step.ep=g;
				step.epp=f;
				step.t=t;
				step.rank=rank;
				step.g_id=-1;
				step.gp_id=gp_id;
				step.gpp_id=gpp_id;
				step.event="S";			      
				sample_steps.push_back(step);
			      }
			    //q_sum+= S_pf_ppg + S_ppf_pg;
			    //S.

			  }

		      //q[g_id][t][e]=q_sum; 

		    }
		  //branches that cross to next time slice  
		  else
		    {
		      //trivial
		      resum=1;
		      if(1)
			{
			  sample_ps.push_back(1);
			  step step;
			  step.e=e;
			  step.ep=-1;
			  step.epp=-1;
			  step.t=t;
			  step.rank=rank;
			  step.g_id=g_id;
			  step.gp_id=-1;
			  step.gpp_id=-1;
			  step.event="0";			
			  sample_steps.push_back(step);      
			}
		      ;//q[g_id][t][e]=q[g_id][t][e];
		    }			  
		}		   
	    }
	  //boundaries for branch e.
	  else
	    {

	      //events within slice rank at time t on branch e 
	      //q[g_id][tpdt][e]=0;
	      //scalar_type q_sum=0;
	      //scalar_type q_sum_nl=0;

	      //non-leaf directed partition		   
	      if (not is_a_leaf)
		for (int i=0;i<N_parts;i++)
		  {	

		    long int gp_id=gp_ids[i];
		    long int gpp_id=gpp_ids[i];	    
		    scalar_type pp=p_part[i];
		    scalar_type qpe=q[gp_id][t][e];
		    scalar_type qppe=q[gpp_id][t][e];
		    scalar_type Sb_pa_ppe= p_Delta_bar*q[gp_id][t][alpha]*qppe*pp;
		    scalar_type Sb_pe_ppa= p_Delta_bar*qpe*q[gpp_id][t][alpha]*pp;
		    //S_bar EVENT
		    resum+=Sb_pa_ppe;
		    if(1)
		      {
			sample_ps.push_back(Sb_pa_ppe);
			step step;
			step.e=-1;
			step.ep=alpha;
			step.epp=e;
			step.t=t;
			step.rank=rank;
			step.g_id=-1;
			step.gp_id=gp_id;
			step.gpp_id=gpp_id;
			step.event="Sb";			      
			sample_steps.push_back(step);
		      }


		    resum+=Sb_pe_ppa;
		    if(1)
		      {
			sample_ps.push_back(Sb_pe_ppa);
			step step;
			step.e=-1;
			step.ep=e;
			step.epp=alpha;
			step.t=t;
			step.rank=rank;
			step.g_id=-1;
			step.gp_id=gp_id;
			step.gpp_id=gpp_id;
			step.event="Sb";		

			sample_steps.push_back(step);	     
		      }
		    //q_sum_nl+= Sb_pa_ppe + Sb_pe_ppa;

		    //q[g_id][tpdt][e]+=p_Delta_bar*(q[gp_id][t][alpha]*q[gpp_id][t][e]+q[gp_id][t][e]*q[gpp_id][t][alpha]);			  
		    //S_bar.

		    scalar_type D=2*p_delta_e*qpe*qppe*pp;
		    resum+=D;
		    if(1)
		      {
			sample_ps.push_back(D);
			step step;
			step.e=-1;
			step.ep=e;
			step.epp=e;
			step.t=t;
			step.rank=rank;
			step.g_id=-1;
			step.gp_id=gp_id;
			step.gpp_id=gpp_id;
			step.event="D";			      

			sample_steps.push_back(step);	     
		      } 
		    //D EVENT
		    //q_sum_nl+= D;

		    //q[g_id][tpdt][e]+=p_delta_e*q[gp_id][t][e]*q[gpp_id][t][e];
		    //D.

		  }

	      scalar_type SLb=p_Delta_bar*Eet*q[g_id][t][alpha];
	      //SL_bar EVENT
	      resum+=SLb;
	      if(1)
		{
		  sample_ps.push_back(SLb);
		  step step;
		  step.e=alpha;
		  step.ep=-1;
		  step.epp=-1;
		  step.t=t;
		  step.rank=rank;
		  step.g_id=g_id;
		  step.gp_id=-1;
		  step.gpp_id=-1;
		  step.event="SLb";			      
		  sample_steps.push_back(step);	     
		}
	      //q_sum_nl+=SLb;

	      //q[g_id][tpdt][e]+=p_Delta_bar*Eet*q[g_id][t][alpha];
	      //SL_bar.

	      //q[g_id][tpdt_nl][e]+=q_sum_nl;

	      scalar_type empty=Get*q[g_id][t][e];
	      //0 EVENT
	      resum+=empty;
	      if(1)
		{
		  sample_ps.push_back(empty);
		  step step;
		  step.e=e;
		  step.ep=-1;
		  step.epp=-1;
		  step.t=t;
		  step.rank=rank;
		  step.g_id=g_id;
		  step.gp_id=-1;
		  step.gpp_id=-1;
		  step.event="0";			   
		  sample_steps.push_back(step);	     
		}
	      //q_sum+=empty;

	      //q[g_id][tpdt][e]=Get*q[g_id][t][e];
	      //0.
		    
	      //q[g_id][tpdt][e]+=q_sum;

	      //events within slice rank at time t on branch e. 	    
	    }
	}
    }

  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################		    
  gp_ids.clear();
  gpp_ids.clear();
  p_part.clear();
  if (S_node)
    {
      rank-=1;      
      t_i=time_slice_times[rank].size()-1;  
      S_node=false;      
    }
  if (set_S_node) 
    {
      S_node=true;
    }
  int step_i=-1;
  scalar_type reresum=0;
  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  scalar_type max_resum=0;
  int max_i=0;
  for (int i=0;i<(int)sample_ps.size();i++)
    {
      if (max_resum<sample_ps[i])
	{
	  max_resum=sample_ps[i];
	  max_i=i;
	}
      reresum+=sample_ps[i];
      if (r*resum<reresum )
	{
	  step_i=i;
	  if (not max_rec) 
	    break;
	}
    }
  if (max_rec)
    step_i=max_i;
  step back_step=sample_steps.at(step_i);
  sample_steps.clear();
  sample_ps.clear();
  /*
    if (back_step.e==-1)
    cout  <<  back_step.event << "\t" << back_step.rank << "\tid_rank:" << id_ranks[back_step.ep] << "\text_sp:" << extant_species[back_step.ep] << "\te:" << back_step.ep<< "\tg_id:" << ale_pointer->set2name(ale_pointer->id_sets[g_id]) << "\t" << back_step.t << " t_i:" << t_i << endl;
    else
    cout  <<  back_step.event << "\t" << back_step.rank << "\tid_rank:" << id_ranks[back_step.e] << "\text_sp:" << extant_species[back_step.e] << "\te:" << back_step.e<< "\tg_id:" << ale_pointer->set2name(ale_pointer->id_sets[g_id]) << "\t" << back_step.t << " t_i:" << t_i << endl;
  */
  stringstream toptmp;
  if (back_step.e==alpha)
    toptmp<<-1;
  else if (id_ranks[back_step.e]==0)
    toptmp<<extant_species[back_step.e];
  else
    toptmp<<id_ranks[back_step.e];
  //cout << branch_length << " +  " << t << " + " << back_step.t << endl;
  scalar_type new_branch_length=-1;//branch_length+t-back_step.t;
  
  size = 0;
  //  if (g_id!=-1) { //We are not at the root bipartition
  
  temp = ale_pointer->id_sets[g_id];
  for (int i = 0; i < ale_pointer->Gamma_size + 1; ++i) {
    // if ( BipartitionTools::testBit ( temp, i) ) {
    if ( temp[i] ) {
      size++;
    }
  }
  //   }

  
  if (ale_pointer->Bip_counts[g_id]>0)
    {
      new_branch_length=max(ale_pointer->Bip_bls[g_id]/ale_pointer->Bip_counts[g_id],(scalar_type)scalar_parameter["min_branch_lenghts"]);
    }
  else
    {
      new_branch_length=max(ale_pointer->Bip_bls[g_id]/ale_pointer->observations,(scalar_type)scalar_parameter["min_branch_lenghts"]);
      
    }

  if (back_step.t==0 and size == 1 and e!=-1)
    {
      register_leaf(e);
      stringstream branch_string;
      if (scalar_parameter["leaf_events"]==1) branch_string<<branch_events;
      branch_string <<":"<<new_branch_length;
      return ale_pointer->set2name(ale_pointer->id_sets[g_id])+branch_string.str();
    }

  if (back_step.event=="D" or back_step.event=="Tb" or back_step.event=="S" or back_step.event=="Sb")
    {

      stringstream transfer_token_stream;
      transfer_token_stream<<"";
      stringstream branch_string;
      if (back_step.event=="S")
	{
	  register_S(e);
	  branch_string<< branch_events
		       <<"."<<id_ranks[e]<<":"<<max(new_branch_length,(scalar_type)0.0); 
	}
      else
	{
	  if (back_step.event=="Tb")
	    {
	      int this_e,this_gid;

	      if (back_step.ep==alpha)
		{
		  this_e=back_step.epp;
		  this_gid=back_step.gpp_id;
		}
	      else
		{
		  this_e=back_step.ep;
		  this_gid=back_step.gp_id;
		}
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
	      tmp<<back_step.rank<<"|"<<t<<"|"<<named_branch.str()<<"|"<<this_gid;
	      register_Ttoken(transfer_token+"|"+tmp.str());
	      // Tto
	      
	      branch_string<< branch_events<<back_step.event<<"@"<<back_step.rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0); 
	    }
	  else if (back_step.event=="Sb")
	    {
	      int this_e;
	      if (back_step.ep==alpha)
		this_e=back_step.epp;
	      else
		this_e=back_step.ep;
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
	      if  (transfer_token!="")
		transfer_token_stream<< transfer_token;
	      else
		transfer_token_stream<< rank<<"|"<<t<<"|"<<named_branch.str()<<"|"<<g_id;
	      branch_string<< branch_events<<"T@"<<rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0); 	    
	    }
	  else
	    {
	      register_D(e);	  
	      stringstream Dtoken_stream;
	      stringstream named_branch;
	      if (e==alpha)
		named_branch<<-1;
	      else if (id_ranks[e]==0)
		named_branch<<extant_species[e];
	      else
		named_branch<<id_ranks[e];

	      Dtoken_stream    << "D|" << rank << "|" <<named_branch.str() << "|"<< g_id;
	      register_Ttoken(Dtoken_stream.str());
	    
	      branch_string<< branch_events<<back_step.event<<"@"<<rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0); 
	    }
	}
      if ( transfer_token_stream.str()=="")
	return "("+
	  sample(S_node,back_step.gp_id, t_i,rank,back_step.ep,0,"",transfer_token,max_rec)
	  +","+
	  sample(S_node,back_step.gpp_id, t_i,rank,back_step.epp,0,"",transfer_token,max_rec)
	  +")"+branch_string.str();
      else
	if(back_step.ep==alpha)
	  return "("+
	    sample(S_node,back_step.gp_id, t_i,rank,back_step.ep,0,"",transfer_token_stream.str(),max_rec)
	    +","+
	    sample(S_node,back_step.gpp_id, t_i,rank,back_step.epp,0,"",transfer_token,max_rec)
	    +")"+branch_string.str();
	else
	  return "("+
	    sample(S_node,back_step.gp_id, t_i,rank,back_step.ep,0,"",transfer_token)
	    +","+
	    sample(S_node,back_step.gpp_id, t_i,rank,back_step.epp,0,"",transfer_token_stream.str(),max_rec)
	    +")"+branch_string.str();
	  
    }
  else if ( back_step.event=="TLb" or back_step.event=="SL" or back_step.event=="SLb" or back_step.event=="0")
    {

      stringstream branch_string;
      stringstream transfer_token_stream;
      transfer_token_stream <<"";
      
      branch_string<< branch_events;
      if (back_step.event!="0")
	{
	  if (back_step.event=="SL")
	    {
	      t_i=time_slice_times[rank].size()-1;
	      register_S(e);	      
	      int f=daughters[e][0];
	      int g=daughters[e][1];
	      if (back_step.e==f)
		register_L(g);
	      else
		register_L(f);
	      branch_string<<"."
			   <<id_ranks[e];
	    }
	  else
	    {
	      if (back_step.event=="TLb")
		{
		  register_Tto(back_step.e);
		  stringstream tmp;

		  stringstream named_branch;
		  if (back_step.e==alpha)
		    named_branch<<-1;
		  else if (id_ranks[back_step.e]==0)
		    named_branch<<extant_species[back_step.e];
		  else
		    named_branch<<id_ranks[back_step.e];

		  tmp<<back_step.rank<<"|"<<t<<"|"<< named_branch.str()<<"|"<<g_id;
		  register_Ttoken(transfer_token+"|"+tmp.str());
		  transfer_token="";
	      
		  branch_string<<""
			       <<"@"<<back_step.rank<<"|"<<named_branch.str();
		}
	      else  if (back_step.event=="SLb")
		{
		  register_L(e);
		  register_Tfrom(e);
		  stringstream named_branch;
		  if (e==alpha)
		    named_branch<<-1;
		  else if (id_ranks[e]==0)
		    named_branch<<extant_species[e];
		  else
		    named_branch<<id_ranks[e];

		  transfer_token_stream<< rank<<"|"<<t<<"|"<<named_branch.str()<<"|"<<g_id;

		  branch_string<<".T"
			       <<"@"<<rank<<"|"<<named_branch.str();
		}
	    }
	}
      if (transfer_token_stream.str()=="")
	{	  
	  return sample(S_node,back_step.g_id, t_i,rank, back_step.e, new_branch_length, branch_string.str(),transfer_token,max_rec);
	}
      else
	{
	  return sample(S_node,back_step.g_id, t_i,rank, back_step.e, new_branch_length, branch_string.str(),transfer_token_stream.str(),max_rec);	
	}
    }
  else
    {
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


