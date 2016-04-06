#include "exODT.h"
using namespace std;
using namespace bpp;

#include <bitset>


static double EPSILON = numeric_limits< double >::min();

//static double EPSILON = 10^-300;

//p(ale) calculates Pi(Gamma) cf. ALEPAPER
scalar_type exODT_model::p(approx_posterior *ale)
{
  ale_pointer=ale;
  //directed partitions and their sizes
  vector <long int>  g_ids;//del-loc
  vector <long int>  g_id_sizes;//del-loc  

  //We sort the directed partitions by size (number of gene tree leaves) to ensure that we calculate things in the proper order (smaller to larger)
  for (map <int, vector <long int > > :: iterator it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	g_ids.push_back((*jt));
	g_id_sizes.push_back((*it).first);
      }
  //root bipartition needs to be handled separately
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
	      p_part.push_back(ale->p_dip(g_id,gp_id,gpp_id));
	  }
      else
	{
	  //root bipartition needs to be handled separately
        map<set<long int>,int> bip_parts;
        for (map <long int,scalar_type> :: iterator it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
	    {
            long int gp_id=(*it).first;
            boost::dynamic_bitset<> gamma =ale->id_sets.at(gp_id);
            boost::dynamic_bitset<> not_gamma = ~gamma;
            not_gamma[0] = 0;
            long int gpp_id = ale->set_ids.at(not_gamma);

            set <long int> parts;
            parts.insert(gp_id);
            parts.insert(gpp_id);
            bip_parts[parts]=1;
            // gamma.clear();
            // not_gamma.clear();
	    }
	  for (map<set<long int>,int> :: iterator kt = bip_parts.begin();kt!=bip_parts.end();kt++)
	    {
	      vector <long int> parts;
            for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) {
                parts.push_back((*sit));
            }
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
		//top of root stem
		{
		  tpdt=t_begin[time_slices[rank][0]];
		  tpdt_rank=rank;
		  tpdt_t_i=0;
		}
	      	     
	      //root
	      scalar_type Delta_t=(tpdt-t)*1;
	      //Delat_bar corresponds to \hat \sigma 
	      scalar_type ni=time_slices[rank].size();
	      scalar_type delta_avg=scalar_parameter["delta_avg"];	      
	      scalar_type tau_avg=scalar_parameter["tau_avg"];	      
	      scalar_type lambda_avg=scalar_parameter["lambda_avg"]; 
	      scalar_type sigma_hat=scalar_parameter["sigma_hat"];
	      scalar_type H_hat=Ee[-1][t];
		    
	      //boundaries for branch alpha virtual branch  

	      //boundary at present
	      if (t==0)
		qvec[g_id+1][rank][t_i][alpha]=0;

	      //boundary between slice rank and rank-1 slice is trivial	
	      ;//qvec[g_id+1][rank][t_i][alpha]=qvec[g_id+1][rank][t_i][alpha];	  

	      //boundaries for branch alpha virtual branch.  
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
			    
			      scalar_type q_sum=0;
			      //qvec[g_id+1][rank][t_i][e]=0;
			    
			      scalar_type SL_fLg=qvec[g_id+1][rank][t_i][f]*Egt;
			      scalar_type SL_Lfg=qvec[g_id+1][rank][t_i][g]*Eft;
			      //SL EVENT, events #3 and #4 in part c of Fig.A1 in http://arxiv.org/abs/1211.4606
			      //qvec[g_id+1][rank][t_i][e]=qvec[g_id+1][rank][t_i][f]*Egt + qvec[g_id+1][rank][t_i][g]*Eft;
			      q_sum+=SL_fLg+SL_Lfg;
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
				    //S EVENT, events #1 and #2 in part c of Fig.A1 in http://arxiv.org/abs/1211.4606
				    //qvec[g_id+1][rank][t_i][e]+=qvec[gp_id+1][rank][t_i][f]*qvec[gpp_id+1][rank][t_i][g] +qvec[gpp_id+1][rank][t_i][f]*qvec[gp_id+1][rank][t_i][g];
				    q_sum+= S_pf_ppg + S_ppf_pg;
				    //S.
				  }
				  //UNDERFLOW ?
				  if (q_sum < EPSILON) {
                    qvec[g_id+1][rank][t_i][e] = EPSILON;
                  }
                  else {
                    qvec[g_id+1][rank][t_i][e] = q_sum;
                  }
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

		  //events within slice rank at time t on alpha virtual branch
		  //scalar_type G_bar=Ge[-1][t];
		  //note that the coalescent approximation in http://arxiv.org/abs/1211.4606 is exp(-(Delta_bar*(n-N)/N+Lambda_bar)*Delta_t );	

		  qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]=0;
		  scalar_type q_sum=0;
		  /* vanishes in the scaling limit  
		  for (int branch_i=0;branch_i<n;branch_i++)			  
		    {
		      int e = time_slices[rank][branch_i];		
		      scalar_type tau_e=vector_parameter["tau"][e];
			  
		      //non-leaf directed partition
		      if (not is_a_leaf)
			for (int i=0;i<N_parts;i++)
			  {	
			    long int gp_id=gp_ids[i];
			    long int gpp_id=gpp_ids[i];	    
			    scalar_type pp=p_part[i];				    
			    scalar_type T_ep_app=p_Ntau_e*qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]*pp;
			    scalar_type T_ap_epp=p_Ntau_e*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]*pp;
			    //T EVENT, events #3 and #4 in part b of Fig.A1 in http://arxiv.org/abs/1211.4606 
			    //qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=p_Ntau_e*(qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]+qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]);
			    q_sum+=T_ep_app+T_ap_epp;
			    //T.

			  }
		    }
		      */

		  //non-leaf directed partition
		  if (not is_a_leaf)
		    for (int i=0;i<N_parts;i++)
		      {	
			long int gp_id=gp_ids[i];
			long int gpp_id=gpp_ids[i];	    
			scalar_type pp=p_part[i];
			
			scalar_type Sb=2*sigma_hat*Delta_t*(qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha])*pp;
			//S_bar EVENT, event #2 in part b of Fig.A1 in http://arxiv.org/abs/1211.4606
			//(note that Delta_bar corresponds to sigma, the Delta_bar,Lambda_bar distinction keeps track of speciaiton (birth) vs extiction (death), 
			// but for the Moran process Delta_bar=Lambda_bar=sigma )
			//qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=p_Delta_bar*(2*qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha]);
			q_sum+=Sb;
			//S_bar.

		      }	    
		

		  for (int branch_i=0;branch_i<n;branch_i++)			  
		    {
		      int e = time_slices[rank][branch_i];		
		      scalar_type tau_e=vector_parameter["tau"][e];
		      scalar_type TLb=tau_e*Delta_t*qvec[g_id+1][rank][t_i][e];
		      //TL_bar EVENT, event #5 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
		      //(note that since Ebar ~ 1, most transfers are expected to involve the TL evenet not the T event,
		      //this should not be confused with the TL event of the Tofigh/Doyon/ODTL models, which here corresponds   
		      // to SL_bar + TL ..)
		      //qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=p_Ntau_e*Ebar*qvec[g_id+1][rank][t_i][e];
		      q_sum+=TLb;
		      //TL_bar.

		    }
		  scalar_type empty=(1+(delta_avg+tau_avg-lambda_avg-ni*sigma_hat-2*sigma_hat*H_hat)*Delta_t)*qvec[g_id+1][rank][t_i][alpha]; 
		  //0 EVENT, event #1 in part b of Fig.A1 in http://arxiv.org/abs/1211.4606
		  //qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=G_bar*qvec[g_id+1][rank][t_i][alpha];
		  q_sum+=empty;
		  //0.

         //UNDERFLOW ?
         qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha]+=q_sum;
         if ( qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha] < EPSILON ) {
            qvec[g_id+1][tpdt_rank][tpdt_t_i][alpha] = EPSILON;
         }

		  //events within slice rank at time t on alpha virtual branch.
		}
	      if(1)
		{
 		  for (int branch_i=0;branch_i<n;branch_i++)
		    {	    
		      int e = time_slices[rank][branch_i];
		      scalar_type Eet=Ee[e][t];	
		      scalar_type delta_e=vector_parameter["delta"][e];
		      scalar_type lambda_e=vector_parameter["lambda"][e];

		      //events within slice rank at time t on branch e 
		      qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=0;
		      scalar_type q_sum=0;

		      //non-leaf directed partition		   
		      if (not is_a_leaf)
			for (int i=0;i<N_parts;i++)
			  {	
			    long int gp_id=gp_ids[i];
			    long int gpp_id=gpp_ids[i];	    
			    scalar_type pp=p_part[i];
			    scalar_type qpe=qvec[gp_id+1][rank][t_i][e];
			    scalar_type qppe=qvec[gpp_id+1][rank][t_i][e];
			    scalar_type Sb_pa_ppe= sigma_hat*Delta_t*qvec[gp_id+1][rank][t_i][alpha]*qppe*pp;
			    scalar_type Sb_pe_ppa= sigma_hat*Delta_t*qpe*qvec[gpp_id+1][rank][t_i][alpha]*pp;

			    //S_bar EVENT, events #3 and #4 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
			    //(The majority of transfer events involve this event.)
			    //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=p_Delta_bar*(qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][e]+qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][alpha]);			  
			    q_sum+= Sb_pa_ppe + Sb_pe_ppa;
			    //S_bar.

			    scalar_type D=2*delta_e*Delta_t*qpe*qppe*pp;
			    //D EVENT, event #2 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
			    //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=p_delta_e*qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][e];
			    q_sum+= D;
			    //D.

			  }

		      scalar_type SLb=sigma_hat*Delta_t*Eet*qvec[g_id+1][rank][t_i][alpha];
		      //SL_bar EVENT, event #5 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
		      //(Transfer events where the donor copy is lost involve this event.)
		      //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=p_Delta_bar*Eet*qvec[g_id+1][rank][t_i][alpha];
		      q_sum+=SLb;
		      //SL_bar.

		      scalar_type empty=(1+(2*delta_e*Eet-sigma_hat*H_hat-delta_e-lambda_e)*Delta_t)*qvec[g_id+1][rank][t_i][e];
		      //0 EVENT, event #1 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
		      //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=Get*qvec[g_id+1][rank][t_i][e];
		      q_sum+=empty;
		      //0.
		   
              qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum;
              //UNDERFLOW?
              if ( qvec[g_id+1][tpdt_rank][tpdt_t_i][e] < EPSILON ) {
                qvec[g_id+1][tpdt_rank][tpdt_t_i][e] = EPSILON;
              }
              
                //if (qvec[g_id+1][tpdt_rank][tpdt_t_i][e]>1) qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=1;
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
  
  scalar_type survive=0;
  scalar_type root_sum=0;
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{
	  scalar_type t=time_slice_times[rank][t_i];			  
	  scalar_type tpdt;
	  if ( t_i < (int)time_slice_times[rank].size()-1 )
	      tpdt=time_slice_times[rank][t_i+1];
	  else if (rank<last_rank-1)
	      tpdt=time_slice_times[rank+1][0];
	  else
	    //top of root stem
	      tpdt=t_begin[time_slices[rank][0]];
	  
	  //root
	  scalar_type Delta_t=(tpdt-t)*1;
	  
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {
	      int e = time_slices[rank][branch_i];
	      //if (rank==last_rank-1 and t_i==(int)time_slice_times[rank].size()-1)//(1-Ee[e][time_slice_times[rank][t_i]])/
	      root_sum+=qvec[0][rank][t_i][e]*Delta_t;
	      survive+=(1-Ee[e][t])*Delta_t;
	      //cout << t<< " " <<rank << " " << branch_i << " " << qvec[0][rank][t_i][e]*Delta_t << " "<< (1-Ee[-1][t])*Delta_t <<endl; 
	    }	     
	  //if (rank==last_rank-1 and t_i==(int)time_slice_times[rank].size()-1)//(1-Ee[-1][time_slice_times[rank][t_i]]);	       
	  root_sum+=qvec[0][rank][t_i][alpha]*Delta_t; 	      
	  survive+=Ee[-1][t]*Delta_t;
	  //cout << t<< " " <<rank << " " << alpha << " " << qvec[0][rank][t_i][alpha]*Delta_t << " "<< (1-Ee[-1][t])*Delta_t <<endl; 

	}
    }

  //del-locs
  g_ids.clear();
  g_id_sizes.clear();

  return root_sum/survive;	
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
  map<int,scalar_type> E_k1,E_k2,E_k3,E_k4;//del-loc
  

  for (int rank=0;rank<last_rank;rank++) 
    for (int tsi=0;tsi<(int)time_slice_times[rank].size();tsi++)
      {
	map<int,map <scalar_type,scalar_type> > y_E;//del-loc
	map<int,map <int,scalar_type> > iy_E;//del-loc

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
	scalar_type sigma_hat=scalar_parameter["sigma_hat"];//vector_parameter["sigma_hat"][rank];//this we usually set to 1
	scalar_type t=t_e;
	scalar_type tpdt=t_b;
	scalar_type h=(tpdt-t)/scalar_parameter["DD"];
	//scalar_type ti=t;
	scalar_type ni=time_slices[rank].size();
	scalar_type delta_avg=scalar_parameter["delta_avg"];	      
	scalar_type tau_avg=scalar_parameter["tau_avg"];	      
	scalar_type lambda_avg=scalar_parameter["lambda_avg"]; 
	scalar_type Delta_bar=Delta_bar;
	scalar_type Lambda_bar=Lambda_bar;
	

	for (int ii=0;ii<scalar_parameter["DD"];ii++)
	  {

	    //intial conditions
	    if (ii==0)
	      {
		if ( t==0)
		  Ee[-1][t]=0; // this is \bar H 
		//trivial else Ee[-1][t]=Ee[-1][t];	      
		//y_E[-1][t]=Ee[-1][t];
		iy_E[-1][ii]=Ee[-1][t];

		//Ee_y[-1]=y_E[-1][t];
		Ee_y[-1]=iy_E[-1][ii];
	      
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
		  }
	      }
	    // RK4: 4th order Runge-Kutta for y'=f(y) 
	    // k1 = f(y[n])	      
	    E_k1[-1]=(-sigma_hat)*Ee_y[-1]*Ee_y[-1]*h+(delta_avg+tau_avg-lambda_avg-sigma_hat*ni)*Ee_y[-1]*h;
	    
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type tau_f=vector_parameter["tau"][f];	       
		E_k1[-1]+=tau_f*(1-Ee_y[f])*h;
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];

		// k1 = f(y[n])	      
		E_k1[e]=delta*Ee_y[e]*Ee_y[e]*h+(-sigma_hat*Ee_y[-1]-delta-lambda)*Ee_y[e]*h+lambda*h;
	      }
	    // k2 = f(y[n]+h/2 k1)

	    //Ee_y[-1]=y_E[-1][ti]+1/2.* E_k1[-1];	      
	    Ee_y[-1]=iy_E[-1][ii]+1/2.* E_k1[-1];	      

	    E_k2[-1]=(-sigma_hat)*Ee_y[-1]*Ee_y[-1]*h+(delta_avg+tau_avg-lambda_avg-sigma_hat*ni)*Ee_y[-1]*h;

	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type tau_f=vector_parameter["tau"][f];	       
		E_k2[-1]+=tau_f*(1-Ee_y[f])*h;

	      }
	  	  
	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		
		// k2 = f(y[n]+h/2 k1)
		//Ee_y[e] =y_E[e][ti]+1/2. * E_k1[e];
		Ee_y[e] =iy_E[e][ii]+1/2. * E_k1[e];
		
		E_k2[e]=delta*Ee_y[e]*Ee_y[e]*h+(-sigma_hat*Ee_y[-1]-delta-lambda)*Ee_y[e]*h+lambda*h;
	      }

	    // k3 = f(y[n]+h/2 k2)
	    //Ee_y[-1]=y_E[-1][ti]+1/2.* E_k2[-1];
	    Ee_y[-1]=iy_E[-1][ii]+1/2.* E_k2[-1];

	    E_k3[-1]=(-sigma_hat)*Ee_y[-1]*Ee_y[-1]*h+(delta_avg+tau_avg-lambda_avg-sigma_hat*ni)*Ee_y[-1]*h;
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type tau_f=vector_parameter["tau"][f];	       
		E_k3[-1]+=tau_f*(1-Ee_y[f])*h;
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
 
		// k3 = f(y[n]+h/2 k2)
		//Ee_y[e] =y_E[e][ti]+1/2. * E_k2[e];
		Ee_y[e] =iy_E[e][ii]+1/2. * E_k2[e];

		E_k3[e]=delta*Ee_y[e]*Ee_y[e]*h+(-sigma_hat*Ee_y[-1]-delta-lambda)*Ee_y[e]*h+lambda*h;
	      }

	    // k4 = f(y[n]+h k3)
	    //Ee_y[-1]=y_E[-1][ti]+1* E_k3[-1];
	    Ee_y[-1]=iy_E[-1][ii]+1* E_k3[-1];

	    E_k4[-1]=(-sigma_hat)*Ee_y[-1]*Ee_y[-1]*h+(delta_avg+tau_avg-lambda_avg-sigma_hat*ni)*Ee_y[-1]*h;
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type tau_f=vector_parameter["tau"][f];	       
		E_k4[-1]+=tau_f*(1-Ee_y[f])*h;
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];

		// k4 = f(y[n]+h k3)
		//Ee_y[e] =y_E[e][ti]+1 * E_k3[e];
		Ee_y[e] =iy_E[e][ii]+1 * E_k3[e];

		E_k4[e]=delta*Ee_y[e]*Ee_y[e]*h+(-sigma_hat*Ee_y[-1]-delta-lambda)*Ee_y[e]*h+lambda*h;
	      }	  
	    if (ii==0)
	      iy_E[-1][ii+1]=Ee[-1][t] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);
	    else
	      iy_E[-1][ii+1]=iy_E[-1][ii] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);

	    

	    if (ii==scalar_parameter["DD"]-1)
	      {
		//Ee[-1][tpdt]=y_E[-1][ti+h];
		Ee[-1][tpdt]=iy_E[-1][ii+1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		// y[n+1] = y[n] + h/6 (k1 + 2 k2 + 2 k3 + k4) 
		if (ii==0)
		  iy_E[e][ii+1]=Ee[e][t] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);
		else
		  iy_E[e][ii+1]=iy_E[e][ii] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);

		if (ii==scalar_parameter["DD"]-1)
		  {
		    //Ee[e][tpdt]=y_E[e][ti+h];
		    Ee[e][tpdt]=iy_E[e][ii+1];
		    //if (e<1) cout << e << " " << t << " " << Ee[e][tpdt] << " " << Ee[-1][tpdt] <<endl;

		  }
	      }
	    //ti=ti+h;	  

	  }      
      }
  //del-locs
  Ee_y.clear();
  E_k1.clear();E_k2.clear();E_k3.clear();E_k4.clear();
}
