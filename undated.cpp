#include "exODT.h"
using namespace std;
using namespace bpp;
static double EPSILON = numeric_limits< double >::min();

void exODT_model::construct_undated(string Sstring)
{
  daughter.clear();
  son.clear();
  name_node.clear();
  node_name.clear();
  node_ids.clear();
  id_nodes.clear();
  
  string_parameter["S_un"]=Sstring;
  S=TreeTemplateTools::parenthesisToTree(string_parameter["S_un"],  (string_parameter["BOOT_STRAP_LABLES"]=="yes")
					 );
  S_root = S->getRootNode();
  vector <Node*> nodes = TreeTemplateTools::getNodes(*S_root);

  for (vector <Node * >::iterator it=nodes.begin();it!=nodes.end();it++ ) (*it)->setDistanceToFather(1);
  
  for (vector <Node * >::iterator it=nodes.begin();it!=nodes.end();it++ )
    if ((*it)->isLeaf())
      {
	name_node[(*it)->getName()]=(*it);
	node_name[(*it)]=(*it)->getName();
      }
    else
      {
	vector<string> leafnames=TreeTemplateTools::getLeavesNames(*(*it));
	sort(leafnames.begin(),leafnames.end());
	stringstream name;
	for (vector <string >::iterator st=leafnames.begin();st!=leafnames.end();st++ )
	  name<<(*st)<<".";
	
	name_node[name.str()]=(*it);
	node_name[(*it)]=name.str();

      }
  // register species
  last_branch=0;
  last_leaf=0;

  set <Node*> saw;
  for (map <string,Node *>::iterator  it=name_node.begin();it!=name_node.end();it++ )
    if ((*it).second->isLeaf())
      {
	Node * node = (*it).second;
	extant_species[last_branch]=node->getName();
	node_ids[node]=last_branch;
	id_nodes[last_branch]=node;
	last_branch++;
	last_leaf++;
	saw.insert(node);
	// a leaf
	daughter[last_branch]=-1;
	// a leaf
	son[last_branch]=-1;
      }
  //ad-hoc postorder
  vector<Node*> next_generation;
  for (map <string,Node *>::iterator  it=name_node.begin();it!=name_node.end();it++ )
    if ((*it).second->isLeaf())
      {
	Node * node = (*it).second;
	next_generation.push_back(node);
      }
  while(next_generation.size())
    {
      vector <Node*> new_generation;
      for (vector<Node*>::iterator  it=next_generation.begin();it!=next_generation.end();it++ )
	{
	  Node * node = (*it);
	  if (node->hasFather() )
	    {
	      Node * father=node->getFather();
	      vector <Node *> sons=father->getSons();
	      Node * sister;
	      if (sons[0]==node) sister=sons[1]; else sister=sons[0];
	      	      
	      if (not node_ids.count(father) and saw.count(sister))
		{
		  node_ids[father]=last_branch;
		  id_nodes[last_branch]=father;		  
		  stringstream name;
		  name << last_branch;
		  father->setBranchProperty("ID",BppString(name.str()));

		  last_branch++;
		  
		  saw.insert(father);
		  new_generation.push_back(father);
		}
	    }
	}
      next_generation.clear();
      for (vector<Node*>::iterator  it=new_generation.begin();it!=new_generation.end();it++ )
	next_generation.push_back((*it));
    }
  string_parameter["S_with_ranks"]=TreeTemplateTools::treeToParenthesis(*S,false,"ID");


  for (map <string,Node *>::iterator  it=name_node.begin();it!=name_node.end();it++ )
    if (not (*it).second->isLeaf())
      {
	Node * node = (*it).second;
	vector <Node *> sons=node->getSons();
	daughter[node_ids[node]]=node_ids[sons[0]];
	son[node_ids[node]]=node_ids[sons[1]];
	//cout << node_ids[node] << " => " << node_ids[sons[0]] << " & " << node_ids[sons[1]] << endl;
	//cout << node_name[node] << " => " << node_name[sons[0]] << " & " << node_name[sons[1]] << endl;  

      }
  branch_counts["Os"].clear();
  branch_counts["Ds"].clear();
  branch_counts["Ts"].clear();
  branch_counts["Tfroms"].clear();
  branch_counts["Ls"].clear();
  branch_counts["count"].clear();
  branch_counts["copies"].clear();
  branch_counts["singleton"].clear();

  for (int e=0;e<last_branch;e++)	
    {
      branch_counts["Os"].push_back(0);
      branch_counts["Ds"].push_back(0);
      branch_counts["Ts"].push_back(0);
      branch_counts["Tfroms"].push_back(0);
      branch_counts["Ls"].push_back(0);
      branch_counts["count"].push_back(0);    
      branch_counts["copies"].push_back(0);
      branch_counts["singleton"].push_back(0);

    }
  T_to_from.clear();
  for (int e=0;e<last_branch;e++)
    {
      vector <scalar_type> tmp;
      T_to_from.push_back(tmp);
      for (int f=0;f<last_branch;f++)
	T_to_from[e].push_back(0);
    }	     
  
  
  last_rank=last_branch;
  set_model_parameter("N",1);

}

void exODT_model::calculate_undatedEs()
{
  uE.clear();
  PD.clear();
  PT.clear();
  PL.clear();
  PS.clear();   
  scalar_type P_T=0;
  for (int f=0;f<last_branch;f++)	  
    P_T+=vector_parameter["tau"][f]/(float)last_branch;
  
  for (int e=0;e<last_branch;e++)
    {
      //cout << e << " is " << node_name[id_nodes[e]] << endl;
      scalar_type P_D=vector_parameter["delta"][e];
      scalar_type P_L=vector_parameter["lambda"][e];
      scalar_type P_S=1;
      scalar_type tmp=P_D+P_T+P_L+P_S;
      P_D/=tmp;
      P_L/=tmp;
      P_S/=tmp;
      PD.push_back(P_D);
      PT.push_back(P_T/tmp);
      PL.push_back(P_L);
      PS.push_back(P_S);
      uE.push_back(0);
    }


  mPTE=0;
  for (int i=0;i<4;i++)
    {
      scalar_type newmPTE=0;
      for (int e=0;e<last_branch;e++)
	{
	  if (e<last_leaf)
	    uE[e]=PL[e]+PD[e]*uE[e]*uE[e]+uE[e]*mPTE;
	  else
	    {
	      int f=daughter[e];
	      int g=son[e];
	      uE[e]=PL[e] + PS[e]*uE[f]*uE[g] + PD[e]*uE[e]*uE[e] + uE[e]*mPTE;
	    }
	  newmPTE+=(PT[e]/(float)last_branch) *uE[e];
	}
      mPTE=newmPTE;
    }
}

scalar_type exODT_model::pun(approx_posterior *ale)
{
  scalar_type survive=0;
  scalar_type root_sum=0;
  uq.clear();mPTuq.clear();//XX
  ale_pointer=ale;

  for (std::map<long int, std::map< scalar_type, std::map<int, scalar_type> > >::iterator it=q.begin();it!=q.end();it++)
    {
      for ( std::map< scalar_type, std::map<int, scalar_type> >::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	(*jt).second.clear();
      (*it).second.clear();
    }      
  q.clear();


  //directed partitions and their sizes
  //vector <long int>  g_ids;
  //vector <long int>  g_id_sizes;
  g_ids.clear();
  g_id_sizes.clear();
  
  for (map <int, vector <long int > > :: iterator it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	g_ids.push_back((*jt));
	g_id_sizes.push_back((*it).first);
      }
  //root bipartition needs to be handled separately
  g_ids.push_back(-1);
  g_id_sizes.push_back(ale->Gamma_size);

  root_i=g_ids.size()-1;
  // gene<->species mapping
  for (int i=0;i<(int)g_ids.size();i++)
    {
      long int g_id=g_ids[i];		
       
      if (g_id_sizes[i]==1)
	{
	  int id = 0;
	  for (auto i=0; i< ale->Gamma_size + 1; ++i)
	    {
	      if ( ale->id_sets[g_id][i] )
		{
		  id=i;
		  break;
		}
	    }
	  
	  string gene_name=ale->id_leaves[ id ];
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
  
  //map <long int, int> g_id2i;
  for (int i=0;i<(int)g_ids.size();i++)
    {
      long int g_id=g_ids[i];	
      g_id2i[g_id]=i;

      if (not ( i<(int)uq.size() ) )
	{
	  vector <scalar_type> tmp;
	  uq.push_back(tmp);
	  mPTuq.push_back(0);       	
	}
      else
	mPTuq[i]=0;
      
      for (int e=0;e<last_branch;e++)
	if (not ( e<(int)uq[i].size() ) )
	  {
	    uq[i].push_back(0);
	  }
	else
	  uq[i][e]=0;
      
    }
  

  for (int iter=0;iter<4;iter++)
    {
      for (int i=0;i<(int)g_ids.size();i++)
	{

	  scalar_type new_mPTuq=0;
	  
	  // directed partition (dip) gamma's id  
	  bool is_a_leaf=false;
	  long int g_id=g_ids[i];	
	  if (g_id_sizes[i]==1)
	    is_a_leaf=true;
	  vector <int> gp_is;
	  vector <long int> gpp_is;
	  vector <scalar_type> p_part;
	  if (g_id!=-1)       
	    for (unordered_map< pair<long int, long int>,scalar_type> :: iterator kt = ale->Dip_counts[g_id].begin(); kt != ale->Dip_counts[g_id].end(); kt++)
	      {	  
		pair<long int, long int> parts = (*kt).first;
		long int gp_id=parts.first;
		long int gpp_id=parts.second;
		int gp_i=g_id2i[parts.first];
		int gpp_i=g_id2i[parts.second];
		gp_is.push_back(gp_i);
		gpp_is.push_back(gpp_i);
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
		  //long int gpp_id=parts[1];

		  int gp_i=g_id2i[parts[0]];
		  int gpp_i=g_id2i[parts[1]];
		  gp_is.push_back(gp_i);
		  gpp_is.push_back(gpp_i);
	      
		  //Here we can create a new ale->Bip_counts[gp_id], in particular for leaves.
		  //We may want to add the leaf entries for Bip_counts when Bip_counts is first created.
		  if (ale->Bip_counts[gp_id]<=scalar_parameter.at("min_bip_count") and not ale->Gamma_size<4)
		    p_part.push_back(0);	      
		  else
		    p_part.push_back(ale->p_bip(gp_id));	      
		}
	      bip_parts.clear();
	    }
	  //######################################################################################################################
	  //#########################################INNNER LOOP##################################################################
	  //######################################################################################################################		    

	  for (int e=0;e<last_branch;e++)
	    {
	      scalar_type uq_sum=0;
	      // S leaf and G leaf 
	      if (e<last_leaf and is_a_leaf and extant_species[e]==gid_sps[g_id])
		{
		  // present
		  uq_sum+=PS[e]*1;
		}
	      // G internal		  
	      if (not is_a_leaf)
		{
		  int N_parts=gp_is.size();		  
		  for (int i=0;i<N_parts;i++)
		    {	
		      int gp_i=gp_is[i];
		      int gpp_i=gpp_is[i];	    
		      scalar_type pp=p_part[i];
		      if (not (e<last_leaf) )
			{
			  int f=daughter[e];
			  int g=son[e];
			  // S event
			  uq_sum+=PS[e]*(uq[gp_i][f]*uq[gpp_i][g] + uq[gp_i][g]*uq[gpp_i][f])*pp;
			}
		      // D event
		      uq_sum+=PD[e]*(uq[gp_i][e]*uq[gpp_i][e]*2)*pp;			      
		      // T event
		      uq_sum+=(uq[gp_i][e]*mPTuq[gpp_i] + uq[gpp_i][e]*mPTuq[gp_i])*pp;			      				      
		    }
		}
	      if (not (e<last_leaf) )
		{
		  int f=daughter[e];
		  int g=son[e];
		  // SL event
		  uq_sum+=PS[e]*(uq[i][f]*uE[g] + uq[i][g]*uE[f]);
		}
	      // DL event
	      uq_sum+=PD[e]*(uq[i][e]*uE[e]*2);			      
	      // TL event
	      uq_sum+=(mPTuq[i]*uE[e] + uq[i][e]*mPTE);			      				      	      
	      if (uq_sum<EPSILON) uq_sum=EPSILON; 
	      uq[i][e]=uq_sum;
	      new_mPTuq +=(PT[e]/(float)last_branch)*uq_sum;
	    }
	  mPTuq[i]=new_mPTuq;
	  
	  //######################################################################################################################
	  //#########################################INNNER LOOP##################################################################
	  //######################################################################################################################		    
	}
      survive=0;
      root_sum=0;
      for (int e=0;e<last_branch;e++)
	{
	  root_sum+=uq[root_i][e];
	  survive+=(1-uE[e]);
	}	     
      
      //cout << root_sum/survive << endl;
    }
  return root_sum/survive;
  
}

string exODT_model::sample_undated()
{

  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  
  scalar_type root_sum=0;
  
  for (int e=0;e<last_branch;e++)
    root_sum+=uq[g_ids.size()-1][e]+EPSILON;
  scalar_type root_resum=0;
  
  for (int e=0;e<last_branch;e++)
    {
      root_resum+=uq[root_i][e]+EPSILON;
      if (r*root_sum<root_resum)
	{
	  register_O(e);
	  return sample_undated(e,root_i,"O")+";";
	}
    }
  return "-!=-";
}

string exODT_model::sample_undated(int e, int i,string last_event,string branch_string)
{

  
  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);

  bool is_a_leaf=false;
  long int g_id=g_ids[i];	
  if (g_id_sizes[i]==1)
    is_a_leaf=true;

  stringstream bl;
  if (ale_pointer->Bip_counts.count(g_id) and ale_pointer->Bip_counts[g_id]>0)
    bl<<max(ale_pointer->Bip_bls[g_id]/ale_pointer->Bip_counts[g_id],(scalar_type)scalar_parameter["min_branch_lenghts"]);
  else
    bl<<max(ale_pointer->Bip_bls[g_id]/ale_pointer->observations,(scalar_type)scalar_parameter["min_branch_lenghts"]); 
  string branch_length=bl.str();

  
  vector <int> gp_is;
  vector <long int> gpp_is;
  vector <scalar_type> p_part;
  if (g_id!=-1)       
    for (unordered_map< pair<long int, long int>,scalar_type> :: iterator kt = ale_pointer->Dip_counts[g_id].begin(); kt != ale_pointer->Dip_counts[g_id].end(); kt++)
      {	  
	pair<long int, long int> parts = (*kt).first;
	long int gp_id=parts.first;
	long int gpp_id=parts.second;
	int gp_i=g_id2i[parts.first];
	int gpp_i=g_id2i[parts.second];
	gp_is.push_back(gp_i);
	gpp_is.push_back(gpp_i);
	if (ale_pointer->Bip_counts[g_id]<=scalar_parameter["min_bip_count"])
	  p_part.push_back(0);
	else
	  p_part.push_back(ale_pointer->p_dip(g_id,gp_id,gpp_id));
      }
  else
    {
      //root bipartition needs to be handled separately
      map<set<long int>,int> bip_parts;
      for (map <long int,scalar_type> :: iterator it = ale_pointer->Bip_counts.begin(); it != ale_pointer->Bip_counts.end(); it++)
	{
	  long int gp_id=(*it).first;
	  boost::dynamic_bitset<> gamma =ale_pointer->id_sets.at(gp_id);
	  boost::dynamic_bitset<> not_gamma = ~gamma;
	  not_gamma[0] = 0;
	  long int gpp_id = ale_pointer->set_ids.at(not_gamma);
	      
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
	  //long int gpp_id=parts[1];

	  int gp_i=g_id2i[parts[0]];
	  int gpp_i=g_id2i[parts[1]];
	  gp_is.push_back(gp_i);
	  gpp_is.push_back(gpp_i);
	      
	  //Here we can create a new ale->Bip_counts[gp_id], in particular for leaves.
	  //We may want to add the leaf entries for Bip_counts when Bip_counts is first created.
	  if (ale_pointer->Bip_counts[gp_id]<=scalar_parameter.at("min_bip_count") and not ale_pointer->Gamma_size<4)
	    p_part.push_back(0);	      
	  else
	    p_part.push_back(ale_pointer->p_bip(gp_id));	      
	}
      bip_parts.clear();
    }

  
  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################		    
  scalar_type uq_sum=0;
  // S leaf and G leaf 
  if (e<last_leaf and is_a_leaf and extant_species[e]==gid_sps[g_id])
    {
      // present
      uq_sum+=PS[e]*1+EPSILON;
    }
  // G internal		  
  if (not is_a_leaf)
    {
      int N_parts=gp_is.size();		  
      for (int i=0;i<N_parts;i++)
	{	
	  int gp_i=gp_is[i];
	  int gpp_i=gpp_is[i];	    
	  scalar_type pp=p_part[i];
	  if (not (e<last_leaf) )
	    {
	      int f=daughter[e];
	      int g=son[e];
	      // S event
	      uq_sum+=PS[e]*uq[gp_i][f]*uq[gpp_i][g]*pp+EPSILON;
	      uq_sum+=PS[e]*uq[gp_i][g]*uq[gpp_i][f]*pp+EPSILON;
	    }
	  // D event
	  uq_sum+=PD[e]*(uq[gp_i][e]*uq[gpp_i][e]*2)*pp+EPSILON;			      
	  // T event
	  for (int f=0;f<last_branch;f++)
	    {
	      uq_sum+=uq[gp_i][e]*(PT[f]/(float)last_branch)*uq[gpp_i][f]*pp+EPSILON;
	      uq_sum+=uq[gpp_i][e]*(PT[f]/(float)last_branch)*uq[gp_i][f]*pp+EPSILON;
	    }
	}
    }
  if (not (e<last_leaf) )
    {
      int f=daughter[e];
      int g=son[e];
      // SL event
      uq_sum+=PS[e]*uq[i][f]*uE[g]+EPSILON;
      uq_sum+=PS[e]*uq[i][g]*uE[f]+EPSILON;
    }
  // DL event
  uq_sum+=PD[e]*(uq[i][e]*uE[e]*2)+EPSILON;			      
  // TL event
  for (int f=0;f<last_branch;f++)
    {
      uq_sum+=(PT[f]/(float)last_branch)*uq[i][f]*uE[e]+EPSILON;
      uq_sum+=(PT[f]/(float)last_branch)*uE[f]*uq[i][e]+EPSILON;			      				      	      
    }
  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################		    

  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################		   
  stringstream estring;
  if (not (e<last_leaf)) estring << e; else estring << extant_species[e]; 
  string estr=estring.str();
    
  scalar_type uq_resum=0;
  // S leaf and G leaf 
  if (e<last_leaf and is_a_leaf and extant_species[e]==gid_sps[g_id])
    {
      // present
      uq_resum+=PS[e]*1+EPSILON;
      if (r*uq_sum<uq_resum)
	{
	  register_leafu(e,last_event);
	  return ale_pointer->set2name(ale_pointer->id_sets[g_id])+branch_string+":"+branch_length;
	}
    }
  // G internal		  
  if (not is_a_leaf)
    {
      int N_parts=gp_is.size();		  
      for (int i=0;i<N_parts;i++)
	{	
	  int gp_i=gp_is[i];
	  int gpp_i=gpp_is[i];	    
	  scalar_type pp=p_part[i];
	  if (not (e<last_leaf) )
	    {
	      int f=daughter[e];
	      int g=son[e];
	      // S event
	      uq_resum+=PS[e]*uq[gp_i][f]*uq[gpp_i][g]*pp+EPSILON;
	      if (r*uq_sum<uq_resum)
		{
		  register_Su(e,last_event);			  
		  return "("+sample_undated(f,gp_i,"S")+","+sample_undated(g,gpp_i,"S")+")."+estr+branch_string+":"+branch_length;
		}		  
	      uq_resum+=PS[e]*uq[gp_i][g]*uq[gpp_i][f]*pp+EPSILON;
	      if (r*uq_sum<uq_resum)
		{
		  register_Su(e,last_event);			  
		  return "("+sample_undated(g,gp_i,"S")+","+sample_undated(f,gpp_i,"S")+")."+estr+branch_string+":"+branch_length;
		}
	    }
	  // D event
	  uq_resum+=PD[e]*(uq[gp_i][e]*uq[gpp_i][e]*2)*pp+EPSILON;
	  if (r*uq_sum<uq_resum)
	    {
	      register_D(e);			 
	      return "("+sample_undated(e,gp_i,"D")+","+sample_undated(e,gpp_i,"D")+").D@"+estr+branch_string+":"+branch_length;
	    }		  
	      
	  // T event
	  for (int f=0;f<last_branch;f++)
	    {
	      stringstream fstring;
	      if (not (f<last_leaf)) fstring << f; else fstring << extant_species[f]; 
	      string fstr=fstring.str();	      
		  
	      uq_resum+=uq[gp_i][e]*(PT[f]/(float)last_branch)*uq[gpp_i][f]*pp+EPSILON;
	      if (r*uq_sum<uq_resum)
		{
		  register_Tfrom(e);
		  register_Tto(f);			  
		  return "("+sample_undated(e,gp_i,"S")+","+sample_undated(f,gpp_i,"T")+").T@"+estr+"->"+fstr+branch_string+":"+branch_length;
		}		  	     
	      uq_resum+=uq[gpp_i][e]*(PT[f]/(float)last_branch)*uq[gp_i][f]*pp+EPSILON;
	      if (r*uq_sum<uq_resum)
		{
		  register_Tfrom(e);
		  register_Tto(f);
		  register_T_to_from(e,f);
		  return "("+sample_undated(e,gpp_i,"S")+","+sample_undated(f,gp_i,"T")+").T@"+estr+"->"+fstr+branch_string+":"+branch_length;
		}		  
		  
	    }
	}
    }
  if (not (e<last_leaf) )
    {
      int f=daughter[e];
      int g=son[e];
      // SL event
      uq_resum+=PS[e]*uq[i][f]*uE[g]+EPSILON;
      if (r*uq_sum<uq_resum)
	{
	  register_Su(e,last_event);
	  register_L(g);			  
	  return sample_undated(f,i,"S","."+estr);
	}		  
      uq_resum+=PS[e]*uq[i][g]*uE[f]+EPSILON;
      if (r*uq_sum<uq_resum)
	{
	  register_Su(e,last_event);
	  register_L(f);			  
	  return sample_undated(g,i,"S","."+estr);
	}		  
    }
  // DL event
  uq_resum+=PD[e]*(uq[i][e]*uE[e]*2)+EPSILON;
  if (r*uq_sum<uq_resum)
    {
      return sample_undated(e,i,"S");
    }		  
  // TL event
  for (int f=0;f<last_branch;f++)
    {
      stringstream fstring;
      if (not (f<last_leaf)) fstring << f; else fstring << extant_species[f]; 
      string fstr=fstring.str();	      
      
      uq_resum+=(PT[f]/(float)last_branch)*uq[i][f]*uE[e]+EPSILON;
      if (r*uq_sum<uq_resum)
	{
	  register_Tfrom(e);
	  register_Tto(f);
	  register_T_to_from(e,f);
	  register_L(e); 
	  return sample_undated(f,i,"T",".T@"+estr+"->"+fstr);
	}		  
      uq_resum+=(PT[f]/(float)last_branch)*uE[f]*uq[i][e]+EPSILON;
      if (r*uq_sum<uq_resum)
	{
	  return sample_undated(e,i,"S");
	}		  
    }
  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################
  cout << "sum error!" << endl;
  return "-!=-";
}

string exODT_model::counts_string_undated()
{
  stringstream out;
  for (int e=0;e<last_branch;e++)	
    {	
      bool isleaf=e<last_leaf;
      stringstream named_branch;
      if (e<last_leaf) named_branch << extant_species[e]; else  named_branch << e; 

      if (not isleaf)
	out<< "S_internal_branch\t"<< named_branch.str() << "\t" 
	   << branch_counts["Ds"][e] << "\t"
	   << branch_counts["Ts"][e] << "\t"
	   << branch_counts["Ls"][e] << "\t"
	   << branch_counts["singleton"][e] << "\t"
	   << branch_counts["copies"][e] << "\n";
      else
	out<< "S_terminal_branch\t"<< named_branch.str() << "\t" 
	   << branch_counts["Ds"][e] << "\t"
	   << branch_counts["Ts"][e] << "\t"
	   << branch_counts["Ls"][e] << "\t"
	   << branch_counts["singleton"][e] << "\t"
	   << branch_counts["copies"][e] << "\n";
	
    }  
  return out.str();
}

void exODT_model::register_Su(int e,string last_event)
{
  MLRec_events["S"]+=1;
  if (e>-1) 
    {
      int f=daughter[e];
      int g=son[e];
      if (last_event=="S" or last_event=="O") branch_counts["singleton"].at(e)+=1;
      branch_counts["copies"].at(e)+=1;
      branch_counts["count"].at(f)+=1;
      branch_counts["count"].at(g)+=1;  
    }
}

void exODT_model::register_leafu(int e,string last_event)
{
  if (e>-1)
    {
      branch_counts["copies"].at(e)+=1;
      if (last_event=="S" or  last_event=="O") branch_counts["singleton"].at(e)+=1;
    }
  //MLRec_events["genes"]+=1;
}

void exODT_model::register_T_to_from(int e,int f)
{
  T_to_from[e][f]+=1;
}

string exODT_model::feSPR(int e, int f)
{
  tree_type * newS=TreeTemplateTools::parenthesisToTree(string_parameter["S_un"],  (string_parameter["BOOT_STRAP_LABLES"]=="yes"));
  Node * newS_root = newS->getRootNode();
  vector <Node*> nodes = TreeTemplateTools::getNodes(*newS_root);

  string e_name=node_name[id_nodes[e]];
  string f_name=node_name[id_nodes[f]];;
  
  Node *e_node, *f_node;

  
  for (vector <Node * >::iterator it=nodes.begin();it!=nodes.end();it++ )
    {
      string name_it;
      if ((*it)->isLeaf())
	{
	  name_it=(*it)->getName();
	}
      else
	{
	  vector<string> leafnames=TreeTemplateTools::getLeavesNames(*(*it));
	  sort(leafnames.begin(),leafnames.end());
	  stringstream name;
	  for (vector <string >::iterator st=leafnames.begin();st!=leafnames.end();st++ )
	    name<<(*st)<<".";
	  
	  name_it=name.str();	  
	}
      if (name_it==e_name) e_node=(*it);
      if (name_it==f_name) f_node=(*it);	
    }

  if (e==f) return string_parameter["S_un"]; 

  bool e_below_f=false;
  Node * node;
  node=e_node;
  while (node->hasFather())
    {
      node=node->getFather();
      if (node==f_node) e_below_f=true; 
    }
  if (e_below_f)
    {
      Node * swap_tmp=e_node;
      e_node=f_node;
      f_node=swap_tmp;
    }
  if (f_node->hasFather() and f_node->getFather()==e_node ) return string_parameter["S_un"]; 
  
  Node * f_father=f_node->getFather();
  vector <Node *> f_sons=f_father->getSons();
  Node * f_sister;
  if (f_sons[0]==f_node) f_sister=f_sons[1]; else f_sister=f_sons[0];
  f_father->removeSon(f_sister);      
  if (f_father->hasFather())
    {
      Node * f_grand_father=f_father->getFather();
      f_grand_father->removeSon(f_father);
      f_grand_father->addSon(f_sister);      
    }
  else
    {
      newS->setRootNode(f_sister);
    }
  if (e_node->hasFather())
    {
      Node * e_father=e_node->getFather();
      e_father->removeSon(e_node);
      e_father->addSon(f_father);
    }
  else
    newS->setRootNode(f_father);

  f_father->addSon(e_node);
  
  for (vector <Node * >::iterator it=nodes.begin();it!=nodes.end();it++ ) (*it)->setDistanceToFather(1);
  
  return TreeTemplateTools::treeToParenthesis(*newS,false,"ID");
}

vector<string> exODT_model::NNIs(int e)
{
  vector<string> NNIs;
  int left_e,right_e,f;
  
  Node * root = id_nodes[e];

  if (root->isLeaf()) return NNIs;
  
  vector <Node *> roots_sons=root->getSons();

  right_e=node_ids[roots_sons[0]];
  left_e=node_ids[roots_sons[1]];    

  if (roots_sons[0]->isLeaf())
    ;
  else
    {
      vector <Node *> right_sons=roots_sons[0]->getSons();
      f=node_ids[right_sons[0]];
      NNIs.push_back(feSPR(left_e,f));
      f=node_ids[right_sons[1]];
      NNIs.push_back(feSPR(left_e,f));
    }
  
  if (roots_sons[1]->isLeaf())
    ;
  else
    {
      vector <Node *> left_sons=roots_sons[1]->getSons();
      f=node_ids[left_sons[0]];
      NNIs.push_back(feSPR(right_e,f));
      f=node_ids[left_sons[1]];
      NNIs.push_back(feSPR(right_e,f));
    }
  return NNIs;
  
}
