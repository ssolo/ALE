#include "ALE.h"
using namespace std;
using namespace bpp;

//## aux. functions ##

//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
template<typename Set> set<Set> powerset(const Set& s, size_t n)
{
    typedef typename Set::const_iterator SetCIt;
    typedef typename set<Set>::const_iterator PowerSetCIt;
    set<Set> res;
    if(n > 0) {
        set<Set> ps = powerset(s, n-1);
        for(PowerSetCIt ss = ps.begin(); ss != ps.end(); ss++)
            for(SetCIt el = s.begin(); el != s.end(); el++) {
                Set subset(*ss);
                subset.insert(*el);
                res.insert(subset);
            }
        res.insert(ps.begin(), ps.end());
    } else
        res.insert(Set());
    return res;
}
//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
template<typename Set> set<Set> powerset(const Set& s)
{
    return powerset(s, s.size());
}

//## approx_posterior class ##
approx_posterior::approx_posterior()
{
  //formal constructor must be followed by load state
  ;
}

approx_posterior::approx_posterior(string tree_string)
{
  constructor_string = tree_string;
  construct(constructor_string);
}

void approx_posterior::construct(string tree_string)
{
  t=new boost::timer();

  last_leafset_id=0;
  observations=0;
  name_seperator="+";
  tree_type * tree = TreeTemplateTools::parenthesisToTree(tree_string);//del-loc  
  Node * root = tree->getRootNode();
  vector <Node * > leaves = TreeTemplateTools::getLeaves(*root);//del-loc
  int id=0;
  for (vector <Node * >::iterator it=leaves.begin();it!=leaves.end();it++ )
    {      
      id++;
      string leaf_name=(*it)->getName();
      leaf_ids[leaf_name]=id;
      Gamma.insert(id);
      id_leaves[id]=leaf_name;
    }
  alpha=0;
  beta=0;

  Gamma_size=Gamma.size();
  //maybe should use boost pow
  //number of partitions of Gamma
  K_Gamma=pow(2.,(int)Gamma_size-1)-1;
  //number of unrooted trees on Gamma_size leaves
  if (Gamma_size==2)
    N_Gamma=1;
  else
    N_Gamma=boost::math::double_factorial<scalar_type>(2*Gamma_size-5);

  //del-locs
  leaves.clear();
  delete tree;
}


void approx_posterior::save_state(string fname)
{
  //constructor_string
  ofstream fout( fname.c_str() );
  //must be first!
  fout<< "#constructor_string" << endl;
  boost::trim(constructor_string);
  fout<< constructor_string << endl;

  fout<< "#observations" <<endl;
  fout<< observations << endl;

  fout<< "#Bip_counts" << endl;
  for (map <long int,scalar_type>::iterator it=Bip_counts.begin();it!=Bip_counts.end();it++)    
    fout<<(*it).first<<"\t"<<(*it).second<< endl;

  fout<< "#Bip_bls" << endl;
  for (map <long int,scalar_type>::iterator it=Bip_bls.begin();it!=Bip_bls.end();it++)    
    fout<<(*it).first<<"\t"<<(*it).second<< endl;

  fout<< "#Dip_counts" << endl;
  for (map <long int,map< set<long int>,scalar_type> >::iterator it=Dip_counts.begin();it!=Dip_counts.end();it++)    
    for (map< set<long int>,scalar_type>::iterator  jt=(*it).second.begin();jt!=(*it).second.end();jt++)    
      {
	fout<<(*it).first<<"\t";
	for (set<long int>::iterator kt=(*jt).first.begin();kt!=(*jt).first.end();kt++)    
	  fout<<(*kt)<<"\t";
	fout<<(*jt).second<<endl;
      }

  fout<< "#last_leafset_id" <<endl;
  fout<< last_leafset_id << endl;

  fout<< "#leaf-id" <<endl;
  for (map<string,int>::iterator it=leaf_ids.begin();it!=leaf_ids.end();it++)    
    fout<<(*it).first<<"\t"<<(*it).second<< endl;

  fout<< "#set-id" <<endl;
  for (map <set <int>,long int>::iterator it=set_ids.begin();it!=set_ids.end();it++)    
    {
      fout << (*it).second;
      fout << "\t:";
      for (set< int>::iterator  jt=(*it).first.begin();jt!=(*it).first.end();jt++)    
	fout << "\t" << (*jt);
      fout << endl;
    }

  fout<< "#END" << endl;

}

void approx_posterior::load_state(string fname)
{
  string tree_string;
  ifstream file_stream (fname.c_str());
  string reading="#nothing";
  if (file_stream.is_open())  //  ########## read state ############
    {
      while (! file_stream.eof())
	{
	 string line;
	 getline (file_stream,line);

	 if (boost::find_first(line, "#"))
	   {
	     boost::trim(line);
	     reading=line;
	   }
	 else if (reading=="#constructor_string")
	   {
	     //cout << reading << endl;
	     boost::trim(line);
	     tree_string=line;
	     constructor_string = tree_string;
	     construct(constructor_string);
	     reading="#nothing";
	   }	 
	 else if (reading=="#observations")
	   {
	     boost::trim(line);	    
	     observations=atof(line.c_str());
	   }
	 else if (reading=="#Bip_counts")
	   {
	     //cout << reading << endl;
	     vector<string> tokens;
	     boost::trim(line);	    
	     boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
	     Bip_counts[atol(tokens[0].c_str())]=atof(tokens[1].c_str());	     
	   }
	 else if (reading=="#Bip_bls")
	   {
	     //cout << reading << endl;
	     vector<string> tokens;
	     boost::trim(line);	    
	     boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
	     Bip_bls[atol(tokens[0].c_str())]=atof(tokens[1].c_str());	     
	   }
	 else if (reading=="#Dip_counts")
	   {
	     //cout << reading << endl;
	     vector<string> tokens;
	     boost::trim(line);	    
	     boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
	     set<long int> parts;
	     parts.insert(atoi(tokens[1].c_str()));
	     parts.insert(atoi(tokens[2].c_str()));
	     Dip_counts[atol(tokens[0].c_str())][parts]=atof(tokens[3].c_str());
	   }
	 else if (reading=="#last_leafset_id")
	   {
	     //cout << reading << endl;
	     boost::trim(line);	    
	     last_leafset_id=atol(line.c_str());	     
	   }
	 else if (reading=="#leaf-id")
	   {
	     //cout << reading << endl;
	     vector<string> tokens;
	     boost::trim(line);	    
	     boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
	     int id=atoi(tokens[1].c_str());
	     string leaf_name=tokens[0];
	     leaf_ids[leaf_name]=id;
	     id_leaves[id]=leaf_name;
	   }
	 else if (reading=="#set-id")
	   {
	     //cout << reading << endl;
	     vector<string> fields;
	     boost::trim(line);	    
	     boost::split(fields,line,boost::is_any_of(":"),boost::token_compress_on);
	     boost::trim(fields[0]);	    
	     long int set_id=atol(fields[0].c_str()); 
	     vector<string> tokens;
	     boost::trim(fields[1]);	    
	     boost::split(tokens,fields[1],boost::is_any_of("\t "),boost::token_compress_on);
	     set <int> id_set;
	     for (vector<string>::iterator it=tokens.begin();it!=tokens.end();it++)
		 id_set.insert(atoi((*it).c_str()));
	     set_ids[id_set]=set_id;
	     id_sets[set_id]=id_set;
	   }
	}
    }
  for (map <long int,set <int> >:: iterator it = id_sets.begin(); it != id_sets.end(); it++)
    {
      set_sizes[(*it).first]=(*it).second.size();      
      size_ordered_bips[(*it).second.size()].push_back((*it).first);
    }
}

string approx_posterior::set2name(set<int> leaf_set)
{
  string name=""; 
  for (set<int>::iterator it=leaf_set.begin();it!=leaf_set.end();it++)
    name+=id_leaves[(*it)]+name_seperator; 
  return name.substr(0,name.size()-1);
}

long int approx_posterior::set2id(set<int> leaf_set)
{	
  long int id=set_ids[leaf_set];
  if (!id)
    {
      last_leafset_id++;
      set_ids[leaf_set]=last_leafset_id;  
      // TMP for debug
      //Dip_levels[leaf_set.size()].push_back(last_leafset_id);
      //id2name[last_leafset_id]=set2name(leaf_set);
      id_sets[last_leafset_id]=leaf_set;
      Bip_bls[last_leafset_id]=0;
      return last_leafset_id;
    }
  else
    {
      return id;
    }
}

scalar_type approx_posterior::Bi(int n2)
{
  int n1=Gamma_size-n2;
  if (n2==1 or n1==1)
    return boost::math::double_factorial<scalar_type>(2*Gamma_size-5);
  n1=max(2,n1);
  n2=max(2,n2);
  return boost::math::double_factorial<scalar_type>(2*n1-3)*boost::math::double_factorial<scalar_type>(2*n2-3);
}

scalar_type approx_posterior::Tri(int n2,int n3)
{
  int n1=Gamma_size-n2-n3;
  n1=max(2,n1);
  n2=max(2,n2);
  n3=max(2,n3);
  return boost::math::double_factorial<scalar_type>(2*n1-3)*boost::math::double_factorial<scalar_type>(2*n2-3)*boost::math::double_factorial<scalar_type>(2*n3-3);
}
scalar_type approx_posterior::binomial(int n,int m)
{
  //maybe worth caching 
  return boost::math::binomial_coefficient<scalar_type>(n,m);
}

scalar_type approx_posterior::trinomial(int n1,int n2, int n3)
{
  //(n,m)!=binomial(n+m,m)
  //(n1,n2,n3)!= (n1+n2,n3)! (n1,n2)! = binomial(n1+n2+n3,n3) binomial(n1+n2,n1)
  // binomial(|Gamma|,i+j) binomial(i+j,j)
  //cf. http://mathworld.wolfram.com/MultinomialCoefficient.html
  return binomial(n1+n2+n3,n3)*binomial(n1+n2,n2);
}

scalar_type approx_posterior::p_bip(set<int> gamma)
{
  if (Gamma_size<4)
    return 1;
  long int g_id=set_ids[gamma];
  return p_bip(g_id);
}


scalar_type approx_posterior::p_dip(set<int> gamma,set<int> gammap,set<int> gammapp)
{
  if (Gamma_size<4)
    return 1;
  long int g_id=set_ids[gamma];  
  long int gp_id=set_ids[gammap];
  long int gpp_id=set_ids[gammapp];
  return p_dip( g_id, gp_id, gpp_id);
}

scalar_type approx_posterior::p_bip(long int g_id)
{
  if (Gamma_size<4)
    return 1;
  scalar_type Bip_count=0;
  if (!g_id)
    {
      //never saw gamma in sample
      Bip_count=0;
    }
  else
    {
      Bip_count=Bip_counts[g_id];
    }
  //if ( gamma.size()==1 or (int)gamma.size()==Gamma_size-1) Bip_count=observations;
  if ( set_sizes[g_id]==1 or set_sizes[g_id]==Gamma_size-1) Bip_count=observations;
  
  //return Bip_count / (observations+alpha) + (alpha/N_Gamma*Bi(gamma.size())) / (observations+alpha);
  return Bip_count / (observations+alpha) + (alpha/N_Gamma*Bi(set_sizes[g_id])) / (observations+alpha);
}
scalar_type approx_posterior::p_dip(long int g_id,long int gp_id,long int gpp_id)
{
  if (Gamma_size<4)
    return 1;
  scalar_type beta_switch=1;
  scalar_type Dip_count=0,Bip_count=0;
  if (!g_id)
    {
      //never saw gamma in sample
      beta_switch=0.;
      Bip_count=0;
      Dip_count=0;
    }
  else
    {
      set <long int> parts;
      parts.insert(gp_id);
      parts.insert(gpp_id);
      Bip_count=Bip_counts[g_id];
      Dip_count=Dip_counts[g_id][parts];      
      if (!gp_id or !gpp_id or Dip_count==0)
	{
	  //never saw gammap-gammapp partition in sample
	  Dip_count=0;
	}
    }
  if (set_sizes[g_id]==1 or set_sizes[g_id]==Gamma_size-1) Bip_count=observations;

  return ( Dip_count + (alpha/N_Gamma*Tri(set_sizes[gp_id],set_sizes[gpp_id])) + beta_switch*beta/(pow(2.,set_sizes[g_id]-1)-1) ) / ( Bip_count + (alpha/N_Gamma*Bi(set_sizes[g_id])) + beta_switch*beta );

}

// an unrooted tree given by its Newick string (which can be rooted)
map <set< int>,scalar_type > approx_posterior::recompose(string G_string)
{
  map <set< int>,scalar_type > return_map;
  map <dedge_type, int> dedges;//del-loc
  map <Node * , vector<Node*> > neighbor;//del-loc

  tree_type * G = TreeTemplateTools::parenthesisToTree(G_string);//del-loc

  if (G->isRooted()) G->unroot();
  vector <Node*> nodes= G -> getNodes(); //del-loc

  //Find all directed edges
  for( vector<Node*>::iterator it=nodes.begin(); it!=nodes.end(); it++)
    {
      Node * from = *it;
      if (from -> hasFather())
	{
	  Node * father = from->getFather();
	  neighbor[from].push_back(father);
	}
      if (! from->isLeaf() )
	{
	  vector<Node*> sons=from->getSons(); //del-loc
	  for( vector<Node*>::iterator it_sons=sons.begin(); it_sons!=sons.end(); it_sons++)
	    neighbor[from].push_back((*it_sons));
	  sons.clear();
	}
      for( vector<Node*>::iterator it_tos=neighbor[from].begin(); it_tos!=neighbor[from].end(); it_tos++)
	{
	  dedge_type dedge;
	  dedge.first = from;
	  dedge.second = *it_tos;
	  dedges[dedge]=0;
	}
    }
  nodes.clear();

  map<dedge_type, set<int> > flat_names; //del-loc
  map<dedge_type, scalar_type > q; //del-loc
  // Name all leaves
  for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
    {
		
      dedge_type dedge=(*it).first;
      Node * from = dedge.first;
      Node * to = dedge.second;
      //visit dedges from a leaf as these can be named  
      if (from->isLeaf())
	{
	  flat_names[dedge].insert(leaf_ids[from->getName()]);	  
	  q[dedge]=1;	
	  return_map[flat_names[dedge]]=q[dedge];
	  //mark named
	  dedges[dedge]=-1;
	  //proceed to dedges in next level -  dedges from cherries can now be named and at least one cherry must exist
	  for( vector<Node*>::iterator it_tos=neighbor[to].begin(); it_tos!=neighbor[to].end(); it_tos++)
	    if ((*it_tos)!=from)
	      {
		dedge_type dedge_out;
		dedge_out.first = to;
		dedge_out.second = *it_tos;
		dedges[dedge_out]+=1;
	      }
	}
    }

  bool edges_left=false;	
  for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
    {
      dedge_type dedge=(*it).first;
      if (dedges[dedge]!=-1)
	edges_left=true;	
    }
  while (edges_left)
    {
      for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
	{
	  dedge_type dedge=(*it).first;			
	  //Process edges that can be named
	  if (dedges[dedge]==2)	
	    {
	      Node * from = dedge.first;
	      Node * to = dedge.second;
	      vector <dedge_type> dedges_in; //del-loc	      
	      for( vector<Node*>::iterator it_tos=neighbor[from].begin(); it_tos!=neighbor[from].end(); it_tos++)
		if (*it_tos!=to)
		  {
		    dedge_type dedge_in;
		    dedge_in.first = *it_tos;
		    dedge_in.second = from;
		    dedges_in.push_back(dedge_in);
		  }

	      
	      set <int> leaf_set_in_1=flat_names[dedges_in[0]];
	      set <int> leaf_set_in_2=flat_names[dedges_in[1]];
	      //flat naming
	      for (set<int>::iterator sit=leaf_set_in_1.begin();sit!=leaf_set_in_1.end();sit++)
		flat_names[dedge].insert((*sit));	
	      for (set<int>::iterator sit=leaf_set_in_2.begin();sit!=leaf_set_in_2.end();sit++)
		flat_names[dedge].insert((*sit));		      	       
	      //flat naming
	      //dip_type dip;
	      //dip.first=set2id(flat_names[dedge]);
	      //dip.second.insert(set2id(leaf_set_in_1));
	      //dip.second.insert(set2id(leaf_set_in_2));
	      
	      q[dedge]=q[dedges_in[0]]*q[dedges_in[1]]*p_dip(flat_names[dedge],leaf_set_in_1,leaf_set_in_2);
	      return_map[flat_names[dedge]]=q[dedge];

	      //mark named
	      dedges[dedge]=-1;
	      //proceed to dedges in next level -  new dedges can now be named 
	      for( vector<Node*>::iterator it_tos=neighbor[to].begin(); it_tos!=neighbor[to].end(); it_tos++)
		if ((*it_tos)!=from)
		  {
		    dedge_type dedge_out;
		    dedge_out.first = to;
		    dedge_out.second = *it_tos;
		    dedges[dedge_out]+=1;
		  }
	      dedges_in.clear();
	    }
	}
      edges_left=false;	
      for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
	{
	  dedge_type dedge=(*it).first;
	  if (dedges[dedge]!=-1)
	    edges_left=true;
	}	
    }

  //del-locs
  dedges.clear();
  for(  map <Node * , vector<Node*> >::iterator it=neighbor.begin(); it!=neighbor.end(); it++)
    (*it).second.clear();
  neighbor.clear();
  delete G;
  for(  map <dedge_type , set<int> >::iterator it=flat_names.begin(); it!=flat_names.end(); it++)
    (*it).second.clear();
  flat_names.clear();
  q.clear();
  return return_map;
}

// an unrooted tree given by its Newick string (which can be rooted)
void approx_posterior::decompose(string G_string, set<int> * bip_ids )
{

  //vector <dip_type > return_dips;
  map <dedge_type, int> dedges;//del-loc
  map <Node * , vector<Node*> > neighbor;//del-loc

  tree_type * G = TreeTemplateTools::parenthesisToTree(G_string);//del-loc

  if (G->isRooted()) G->unroot();
  vector <Node*> nodes= G -> getNodes(); //del-loc

  //Find all directed edges
  for( vector<Node*>::iterator it=nodes.begin(); it!=nodes.end(); it++)
    {
      Node * from = *it;
      if (from -> hasFather())
	{
	  Node * father = from->getFather();
	  neighbor[from].push_back(father);
	}
      if (! from->isLeaf() )
	{
	  vector<Node*> sons=from->getSons(); //del-loc
	  for( vector<Node*>::iterator it_sons=sons.begin(); it_sons!=sons.end(); it_sons++)
	    neighbor[from].push_back((*it_sons));
	  sons.clear();
	}
      for( vector<Node*>::iterator it_tos=neighbor[from].begin(); it_tos!=neighbor[from].end(); it_tos++)
	{
	  dedge_type dedge;
	  dedge.first = from;
	  dedge.second = *it_tos;
	  dedges[dedge]=0;
	}
    }
  nodes.clear();

  map<dedge_type, set<int> > flat_names; //del-loc
  // Name all leaves
  for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
    {
		
      dedge_type dedge=(*it).first;
      Node * from = dedge.first;
      Node * to = dedge.second;
      //visit dedges from a leaf as these can be named  
      if (from->isLeaf())
	{
	  flat_names[dedge].insert(leaf_ids[from->getName()]);
	  //bl - hack
	  long int g_id=set2id(flat_names[dedge]);
	  if (from->hasDistanceToFather())
	    Bip_bls[g_id]+=from->getDistanceToFather();		  
	  else
	    Bip_bls[g_id]+=0;		  
	  //mark named
	  dedges[dedge]=-1;
	  //proceed to dedges in next level -  dedges from cherries can now be named and at least one cherry must exist
	  for( vector<Node*>::iterator it_tos=neighbor[to].begin(); it_tos!=neighbor[to].end(); it_tos++)
	    if ((*it_tos)!=from)
	      {
		dedge_type dedge_out;
		dedge_out.first = to;
		dedge_out.second = *it_tos;
		dedges[dedge_out]+=1;
	      }
	}
    }

  bool edges_left=false;	
  for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
    {
      dedge_type dedge=(*it).first;
      if (dedges[dedge]!=-1)
	edges_left=true;	
    }

  if (G -> getLeaves().size()==2)
    {
      Bip_counts[(long int) 1]+=1;
    }
  else
  while (edges_left)
    {
      for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
	{

	  dedge_type dedge=(*it).first;			
	  //Process edges that can be named
	  if (dedges[dedge]==2)	
	    {
	      Node * from = dedge.first;
	      Node * to = dedge.second;
	      vector <dedge_type> dedges_in; //del-loc	      
	      for( vector<Node*>::iterator it_tos=neighbor[from].begin(); it_tos!=neighbor[from].end(); it_tos++)
		if (*it_tos!=to)
		  {
		    dedge_type dedge_in;
		    dedge_in.first = *it_tos;
		    dedge_in.second = from;
		    dedges_in.push_back(dedge_in);
		  }
	      
	      set <int> leaf_set_in_1=flat_names[dedges_in[0]];
	      set <int> leaf_set_in_2=flat_names[dedges_in[1]];
	      //flat naming
	      for (set<int>::iterator sit=leaf_set_in_1.begin();sit!=leaf_set_in_1.end();sit++)
		flat_names[dedge].insert((*sit));	
	      for (set<int>::iterator sit=leaf_set_in_2.begin();sit!=leaf_set_in_2.end();sit++)
		flat_names[dedge].insert((*sit));		      	       
	      //flat naming
	      
	      long int g_id=set2id(flat_names[dedge]);
	      set <long int> parts;
	      parts.insert(set2id(leaf_set_in_1));
	      parts.insert(set2id(leaf_set_in_2));

	      //bl - hack
	      if (from->hasFather() and from->getFather()==to)
		{
		  if (from->hasDistanceToFather())
		    Bip_bls[g_id]+=from->getDistanceToFather();		  
		  else
		    Bip_bls[g_id]+=0;		  
		}
	      else if (to->hasFather() and to->getFather()==from)
		{
		  if (to->hasDistanceToFather())
		    Bip_bls[g_id]+=to->getDistanceToFather();
		  else
		    Bip_bls[g_id]+=0;
		}
	      else
		{
		  cout << "impossible" <<endl;
		}
	      //bl - hack

	      //INTEGRATED COUNTING 
	      Dip_counts[g_id][parts]+=1;
	      Bip_counts[g_id]+=1;
	      //bipartion naming
	      if (bip_ids!=NULL) bip_ids->insert(g_id);
	      
	      //dip.first=g_id;
	      //dip.second=parts;
	      //return_dips.push_back(dip);

	      //mark named
	      dedges[dedge]=-1;
	      //proceed to dedges in next level -  new dedges can now be named 
	      for( vector<Node*>::iterator it_tos=neighbor[to].begin(); it_tos!=neighbor[to].end(); it_tos++)
		if ((*it_tos)!=from)
		  {
		    dedge_type dedge_out;
		    dedge_out.first = to;
		    dedge_out.second = *it_tos;
		    dedges[dedge_out]+=1;
		  }
	      dedges_in.clear();

	    }

	}
      edges_left=false;	
      for( map<dedge_type, int>::iterator it=dedges.begin(); it!=dedges.end(); it++)
	{
	  dedge_type dedge=(*it).first;
	  if (dedges[dedge]!=-1)
	    edges_left=true;
	}	
    }

  //del-locs
  dedges.clear();
  for(  map <Node * , vector<Node*> >::iterator it=neighbor.begin(); it!=neighbor.end(); it++)
    (*it).second.clear();
  neighbor.clear();
  delete G;
  for(  map <dedge_type , set<int> >::iterator it=flat_names.begin(); it!=flat_names.end(); it++)
    (*it).second.clear();
  flat_names.clear();  
  //return return_dips;
}

void approx_posterior::observation(vector<string> trees, bool count_topologies)
{
  for (vector<string>::iterator it=trees.begin();it!=trees.end();it++)
    {
      //cout << (*it) << endl;
      if (count_topologies)
	{
	  set <int> bip_ids;//del-loc      
	  decompose(*it,&bip_ids);//del-loc      
	  string bip_string="|";
	  for (set <int>::iterator st=bip_ids.begin();st!=bip_ids.end();st++)
	    bip_string+=(*st)+"|";
	  if (bipstring_trees.count(bip_string)==0)
	    {
	      bipstring_trees[bip_string]=*it;
	      tree_bipstrings[*it]=bip_string;
	    }
	  tree_counts[bipstring_trees[bip_string]]+=1;
	  bip_ids.clear();
	}
      else decompose(*it);//del-loc            
      observations+=1;
    }
  //cout << "obsdone." << endl;

  
  set_sizes.clear();
  for (map <int, vector <long int > > :: iterator it = size_ordered_bips.begin(); it != size_ordered_bips.end(); it++)
    (*it).second.clear();
  size_ordered_bips.clear();

  for (map <long int,set <int> >:: iterator it = id_sets.begin(); it != id_sets.end(); it++)
    {
      set_sizes[(*it).first]=(*it).second.size();      
      size_ordered_bips[(*it).second.size()].push_back((*it).first);
    }
}


// of an unrooted tree given by its Newick string (which can be rooted)
scalar_type approx_posterior::p(string tree_string)
{
  scalar_type p=0;
  map <set<int>,scalar_type> rec_map=recompose( tree_string);
  for (map <set<int>,scalar_type>::iterator it=rec_map.begin();it!=rec_map.end();it++) 
    {
      p=(*it).second;
      set <int> gamma=(*it).first;
      set <int> not_gamma;
      for (set<int>::iterator st=Gamma.begin();st!=Gamma.end();st++)
	if (gamma.count(*st)==0)
	  not_gamma.insert(*st);
      p*=rec_map[not_gamma]*p_bip(gamma);
      break;
    }
  return p;
}

pair<string,scalar_type> approx_posterior::mpp_tree()
{
  map<long int, scalar_type > qmpp; //del-loc
  for (map <int, vector <long int > > :: iterator it = size_ordered_bips.begin(); it != size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	long int g_id=(*jt);
	// leaves
	if ((*it).first==1)
	  qmpp[g_id]=1;
	else
	  {
	    scalar_type max_cp=0;
	    for (map< set<long int>,scalar_type> :: iterator kt = Dip_counts[g_id].begin(); kt != Dip_counts[g_id].end(); kt++)
	      {	  
		vector <long int> parts;
		for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
		long int gp_id=parts[0];
		long int gpp_id=parts[1];	    
		scalar_type cp=p_dip(g_id,gp_id,gpp_id)*qmpp[gp_id]*qmpp[gpp_id];	
		if (cp>max_cp) max_cp=cp;
	      }
	    qmpp[g_id]=max_cp;
	  }
      }
  scalar_type max_pp=0,sum_pp=0;
  long int max_bip=-1,max_not_bip=-1;
  //we look at everything twice..
  for (map <long int,scalar_type> :: iterator it = Bip_counts.begin(); it != Bip_counts.end(); it++)
    {
      long int g_id=(*it).first;
      set <int> gamma=id_sets[g_id];
      set <int> not_gamma;
      for (set<int>::iterator st=Gamma.begin();st!=Gamma.end();st++)
	if (gamma.count(*st)==0)
	  not_gamma.insert(*st);
      long int not_g_id = set_ids[not_gamma];
      scalar_type pp=qmpp[g_id]*qmpp[not_g_id]*p_bip(g_id);
      sum_pp+=pp;
      if (max_pp<pp) {max_pp=pp; max_bip=g_id; max_not_bip=not_g_id;}
    }    
  stringstream bs;
  //we looked at everything twice..
  bs<<max_pp/sum_pp*2<<":"<<min(Bip_bls[max_bip]/Bip_counts[max_bip],(scalar_type)0.99);

  string max_tree="("+mpp_backtrack(max_bip,&qmpp)+","+mpp_backtrack(max_not_bip,&qmpp)+")"+bs.str()+";\n";

  //cout << max_tree << endl;
  pair <string,scalar_type> return_pair;
  return_pair.first=max_tree;
  return_pair.second=max_pp;
  return return_pair;
}

 string approx_posterior::mpp_backtrack(long int g_id, map<long int, scalar_type > * qmpp)
{
  //leaf
  if (set_sizes[g_id]==1) 
    {
      stringstream bs;
      bs<<Bip_bls[g_id]/observations;
      return id_leaves[(*id_sets[g_id].begin())]+":"+bs.str(); 
    }

  scalar_type max_cp=0,sum_cp=0;
  long int max_gp_id=-1;
  long int max_gpp_id=-1;
  for (map< set<long int>,scalar_type> :: iterator kt = Dip_counts[g_id].begin(); kt != Dip_counts[g_id].end(); kt++)
    {	  
      vector <long int> parts;
      for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
      long int gp_id=parts[0];
      long int gpp_id=parts[1];	    
      scalar_type cp=p_dip(g_id,gp_id,gpp_id)*(*qmpp)[gp_id]*(*qmpp)[gpp_id];		
      sum_cp+=cp;
      if (cp>max_cp) {max_cp=cp; max_gp_id=gp_id; max_gpp_id=gpp_id;}
    }
  stringstream bs;
  bs<<max_cp/sum_cp<<":"<<Bip_bls[g_id]/Bip_counts[g_id];
  return "("+mpp_backtrack(max_gp_id,qmpp)+","+mpp_backtrack(max_gpp_id,qmpp)+")"+bs.str();
}
 
//random tree from unrooted posterior
string approx_posterior::random_tree()
{
  // we start at an observed partition
  
  set<int> gamma;
  scalar_type sum=0;  
  for (map <long int,scalar_type>::iterator it=Bip_counts.begin();it!=Bip_counts.end();it++)    
    {
      sum+=(*it).second;
    }
  scalar_type rnd=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);    
  scalar_type re_sum=0;
  long int g_id;
  for (map <long int,scalar_type>::iterator it=Bip_counts.begin();it!=Bip_counts.end();it++)    
    {
      re_sum+=(*it).second;
      g_id=(*it).first;
      if (re_sum>sum*rnd)
	break;
    }
  gamma = id_sets[g_id];
  set <int> not_gamma;
  for (int i=1;i<Gamma_size+1;i++) if (!gamma.count(i)) not_gamma.insert(i);
  return "("+random_split(gamma)+":1,"+random_split(not_gamma)+":1);\n";
}

string approx_posterior::random_split(set <int> gamma)
{
  // if gamma contains only a leaf we return its name
  vector <int> gamma_v;
  // std::set-s are ordered and SHOULD have random acces, but don't, hence this 
  for (set<int>::iterator sit=gamma.begin();sit!=gamma.end();sit++) gamma_v.push_back((*sit));
  if (gamma.size()==1)
    return id_leaves[ gamma_v[0]];
  scalar_type p_sum=0;  
  //rnd for choosing directed partition 
  scalar_type rnd=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);  
  long int gp_id,gpp_id,g_id;
  scalar_type Bip_count,beta_switch=1;
  g_id=set_ids[gamma];
  set <int> gammap,gammapp;
  if (!g_id)
    {
      //never saw gamma in sample
      beta_switch=0.;
      Bip_count=0;
    }
 
  for (int gp_size=1;gp_size<=(int)gamma.size()/2;gp_size++)
    {
      int saw=0;
      //see if a directed partition is the one we choose
      if (g_id) 
	for (map< set<long int> ,scalar_type> ::iterator dit=Dip_counts[g_id].begin();dit!=Dip_counts[g_id].end();dit++)    
	  { 
	    vector<long int> parts_v;
	    // std::set-s are ordered and SHOULD have random acces, but don't hence this 
	    for (set<long int>::iterator sit=(*dit).first.begin();sit!=(*dit).first.end();sit++) parts_v.push_back((*sit));
	    gp_id=parts_v[0];
	    gpp_id=parts_v[1];
	    parts_v.clear();
	    int this_size=min((int)id_sets[gp_id].size(),(int)id_sets[gpp_id].size());
	    if ( this_size==gp_size )
	      {
		p_sum+=p_dip(gamma,id_sets[gp_id],id_sets[gpp_id]);
		//see if this directed partition is the one we choose
		if (rnd<p_sum)	
		  p_sum=-1;
		saw+=1;
	      }
	    if (p_sum<0)
	      {
		gammap=id_sets[gp_id];
		gammapp=id_sets[gpp_id];
		break;
	      }
	  }      
      if (p_sum<0)
      	break;
      //sum the prob.s of all unobserved bipartitons
      Bip_count=Bip_counts[g_id];
      if (gamma.size()==1 or (int)gamma.size()==Gamma_size-1) Bip_count=observations;
      int nbip=binomial(gamma.size(),gp_size);
      if ((int)gamma.size()-gp_size==gp_size) nbip/=2;
      p_sum+=  ( 0 + (alpha/N_Gamma*Tri(gp_size,gamma.size()-gp_size)) + beta_switch*beta/(pow(2.,(int)gamma.size()-1)-1) ) / ( Bip_count + (alpha/N_Gamma*Bi(gamma.size())) + beta_switch*beta )*(nbip-saw);
      
      //see if an unsampled directed partition is the one we choose
      
      if (rnd<p_sum)	
	p_sum=-1;
      if (p_sum<0)
	{
	  //pick one of these at random
	  bool stop=false;
	  while (!stop)
	    {
	      gammap.clear();gammapp.clear();
	      for (int i=0;i<gp_size;i++)
		gammap.insert( gamma_v[(RandomTools::giveIntRandomNumberBetweenZeroAndEntry(gamma.size()))] );       
	      for (int i=0;i<(int)gamma.size();i++)  if (!gammap.count(gamma_v[i])) gammapp.insert(gamma_v[i]);	    
	      gp_id=set_ids[gammap];
	      gpp_id=set_ids[gammapp];
	      set <long int> parts;
	      parts.insert(gp_id);
	      parts.insert(gpp_id);
	      if (Dip_counts[g_id][parts]==0) stop=true;
	    }
	  break;
	}
      
    }
  return "("+random_split(gammap)+":1,"+random_split(gammapp)+":1)";
      
}


vector<string> approx_posterior::all_trees(set <int> gamma)
{
  vector< string > all_trees_g;
  if (gamma.size()==1)
    {
      all_trees_g.push_back(id_leaves[ *(gamma.begin()) ]);
    }
  else
    {
      set< set<int> > P_gamma=powerset< set<int> >(gamma);//del-loc  
      for (set<set<int> >::iterator st=P_gamma.begin();st!=P_gamma.end();st++)
	{
	  if (gamma.size()>(*st).size() and (*st).size()>0 and (*st).count(*(gamma.begin()))==1)
	    {

	      set <int> not_st;//del-loc
	      
	      vector< string > all_trees_gp=all_trees(*st);//del-loc
	      
	      for (set<int>::iterator nst=gamma.begin();nst!=gamma.end();nst++)
		if ((*st).count(*nst)==0)
		  not_st.insert(*nst);
	      
	      vector< string > all_trees_gpp=all_trees(not_st);//del-loc
	      
	      for (vector<string>::iterator lt=all_trees_gp.begin();lt!=all_trees_gp.end();lt++)
		for (vector<string>::iterator rt=all_trees_gpp.begin();rt!=all_trees_gpp.end();rt++)
		  {
		    all_trees_g.push_back("("+(*lt)+","+(*rt)+")");	
		  }
	      not_st.clear();
	      all_trees_gp.clear();
	      all_trees_gpp.clear();
	    }
	}

      P_gamma.clear();      
    }
  if (Gamma==gamma)
    for (vector<string>::iterator it=all_trees_g.begin();it!=all_trees_g.end();it++)
      {
	(*it)+=";\n";
      }
  return all_trees_g;
}


scalar_type approx_posterior::count_all_trees(set <int> gamma)
{
  scalar_type count_trees_g=0;
  if (gamma.size()==1)
    {
      count_trees_g=1;
    }
  else
    {
      set< set<int> > P_gamma=powerset< set<int> >(gamma);//del-loc  
      for (set<set<int> >::iterator st=P_gamma.begin();st!=P_gamma.end();st++)
	{
	  if (gamma.size()>(*st).size() and (*st).size()>0 and (*st).count(*(gamma.begin()))==1)
	    {

	      set <int> not_st;//del-loc
	      
	      scalar_type count_trees_gp=count_all_trees(*st);//del-loc
	      
	      for (set<int>::iterator nst=gamma.begin();nst!=gamma.end();nst++)
		if ((*st).count(*nst)==0)
		  not_st.insert(*nst);
	      
	      scalar_type count_trees_gpp=count_all_trees(not_st);//del-loc

	      count_trees_g+=count_trees_gp*count_trees_gpp;
	      not_st.clear();
	    }
	}

      P_gamma.clear();      
    }
  return count_trees_g;
}

scalar_type approx_posterior::count_trees()
{
  scalar_type count_trees_g=0;
  
  map<long int,scalar_type> g_id_count;//del-loc

  for (map <int, vector <long int > > :: iterator it = size_ordered_bips.begin(); it != size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	long int g_id=(*jt);
	// leaves
	if ((*it).first==1)
	   g_id_count[g_id]=1;
	else
	  {
	    g_id_count[g_id]=0;
	    for (map< set<long int>,scalar_type> :: iterator kt = Dip_counts[g_id].begin(); kt != Dip_counts[g_id].end(); kt++)
	      {	  
		vector <long int> parts;
		for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
		long int gp_id=parts[0];
		long int gpp_id=parts[1];	    
		g_id_count[g_id]+=g_id_count[gp_id]*g_id_count[gpp_id];
	      }
	  }
      }

  for (map <long int,scalar_type> :: iterator it = Bip_counts.begin(); it != Bip_counts.end(); it++)
    {
      long int g_id=(*it).first;
      set <int> gamma=id_sets[g_id];
      set <int> not_gamma;
      for (set<int>::iterator st=Gamma.begin();st!=Gamma.end();st++)
	if (gamma.count(*st)==0)
	  not_gamma.insert(*st);
      if ( gamma.size()>not_gamma.size())
	count_trees_g+=g_id_count[set_ids[gamma]]*g_id_count[set_ids[not_gamma]];//count_trees(set_ids[gamma])*count_trees(set_ids[not_gamma]);
      else if (gamma.size()==not_gamma.size())
	count_trees_g+=g_id_count[set_ids[gamma]]*g_id_count[set_ids[not_gamma]]/2.0;//count_trees(set_ids[gamma])*count_trees(set_ids[not_gamma])/2.0;
      //cout << count_trees(gamma) << " " << set2name(gamma) << " " << count_trees(not_gamma) << " " << set2name(not_gamma) <<endl; 
    }
  g_id_count.clear();
  return count_trees_g;
}
scalar_type approx_posterior::count_trees(long int g_id)
{
  scalar_type count_trees_g=0;
  //std::map <std::set <int>,long int>  set_ids;//del-loc  
  //std::map< long int, std::set <int> > id_sets;//del-loc
  //long int g_id=set_ids[gamma];
  set<int> gamma=id_sets[g_id];
  int gamma_size=gamma.size();
  if (gamma_size==1)
    {
      count_trees_g=1;
    }
  else
    {
      set< long int > P_gamma;//=powerset< set<int> >(gamma);//del-loc  
      for (map< set<long int>,scalar_type> :: iterator kt = Dip_counts[g_id].begin(); kt != Dip_counts[g_id].end(); kt++)
	{
	  vector <long int> parts;
	  for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
	  long int gp_id=parts[0];
	  long int gpp_id=parts[1];	    
	  P_gamma.insert(gp_id);
	  P_gamma.insert(gpp_id);
	}
      for (set<long int>::iterator st=P_gamma.begin();st!=P_gamma.end();st++)
	{
	  set<int> gammap=id_sets[(*st)];
	  if (gamma.size()>gammap.size() and gammap.size()>0 and gammap.count(*(gamma.begin()))==1)
	    {

	      set <int> not_gammap;//del-loc    
	      scalar_type count_trees_gp=count_trees((*st));//del-loc	      
	      for (set<int>::iterator nst=gamma.begin();nst!=gamma.end();nst++)
		if (gammap.count(*nst)==0)
		  not_gammap.insert(*nst);	      
	      scalar_type count_trees_gpp=count_trees(set_ids[not_gammap]);//del-loc
	      not_gammap.clear();

	      count_trees_g+=count_trees_gp*count_trees_gpp;
	    }
	  gammap.clear();
	}
      gamma.clear();
      P_gamma.clear();      
    }
  return count_trees_g;
}



