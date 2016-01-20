#include "ALE.h"
using namespace std;
using namespace bpp;

int prune(string * tree,vector<string> keep_list,vector<string> * keep_leaves)
{

  map <string,int> keep_map;
  for (vector<string>::iterator it=keep_list.begin();it!=keep_list.end();it++) keep_map[(*it)]=1;    
  tree_type * T=TreeTemplateTools::parenthesisToTree(*tree,false,"ID");
  vector <Node *> leaves=T->getLeaves();
  int N_leaves=leaves.size();
  for (vector <Node *>::iterator it=leaves.begin();it!=leaves.end();it++)
    {
      string name=(*it)->getName();
      boost::trim(name);
      vector<string> tokens;
      boost::split(tokens,name,boost::is_any_of("_"),boost::token_compress_on);
      string sp=tokens[0];
      if (keep_map.count(sp)==1) (*keep_leaves).push_back(name);
      if (not keep_map.count(sp)==1 )
	{
	  if (N_leaves>3) TreeTemplateTools::dropLeaf(*T,name);
	  N_leaves--;
	}
    }    
  if (N_leaves>3)  *tree=TreeTemplateTools::treeToParenthesis(*T,false,"ID");
  return N_leaves;
}

int cout_small(vector <string> keep_leaves,string ale_name)
{
  ofstream fout(ale_name.c_str());
  cout << "#wrting "<< ale_name << endl;
  if (keep_leaves.size()==1)
    {
      fout << "#constructor_string"<< endl;
      fout << keep_leaves[0]<< endl;
      fout << "#observations"<< endl;
      fout << "1"<< endl;
      fout << "#Bip_counts"<< endl;
      fout << "#Bip_bls"<< endl;
      fout << "1 1"<< endl;
      fout << "#Dip_counts"<< endl;
      fout << "#last_leafset_id"<< endl;
      fout << "1"<< endl;
      fout << "#leaf-id"<< endl;
      fout << keep_leaves[0]+" 1"<< endl;
      fout << "#set-id"<< endl;
      fout << "1 : 1"<< endl;
      fout << "#END"<< endl;
    }
  else if (keep_leaves.size()==2)
    {
      fout << "#constructor_string"<< endl;
      fout << keep_leaves[0]+","+keep_leaves[1]<< endl;
      fout << "#observations"<< endl;
      fout << "1"<< endl;
      fout << "#Bip_counts" << endl;
      fout << "1 1"<< endl;
      fout << "2 1"<< endl;
      fout << "#Bip_bls"<< endl;
      fout << "1 1"<< endl;
      fout << "2 1"<< endl;
      fout << "#Dip_counts"<< endl;
      fout << "#last_leafset_id"<< endl;
      fout << "2"<< endl;
      fout << "#leaf-id"<< endl;
      fout << keep_leaves[0]+" 1"<< endl;
      fout << keep_leaves[1]+" 2"<< endl;
      fout << "#set-id"<< endl;
      fout << "1 : 1"<< endl;
      fout << "2 : 2"<< endl;
      fout << "#END"<< endl;
    }
  else if (keep_leaves.size()==3)
    {
      fout << "#constructor_string"<< endl;
      fout << keep_leaves[0]+","+keep_leaves[1]+","+keep_leaves[2]<< endl;
      fout << "#observations"<< endl;
      fout << "1"<< endl;
      fout << "#Bip_counts" << endl;
      fout << "4 1"<< endl;
      fout << "5 1"<< endl;
      fout << "6 1"<< endl;
      fout << "#Bip_bls"<< endl;
      fout << "1 1"<< endl;
      fout << "2 1"<< endl;
      fout << "3 1"<< endl;    
      fout << "4 1"<< endl;
      fout << "5 1"<< endl;
      fout << "6 1"<< endl;
      fout << "#Dip_counts"<< endl;
      fout << "4       2       3       1"<< endl;
      fout << "5       1       3       1"<< endl;
      fout << "6       1       2       1"<< endl;
      fout << "#last_leafset_id"<< endl;
      fout << "6"<< endl;
      fout << "#leaf-id"<< endl;
      fout << keep_leaves[0]+" 1"<< endl;
      fout << keep_leaves[1]+" 2"<< endl;
      fout << keep_leaves[2]+" 3"<< endl;
      fout << "#set-id"<< endl;
      fout << "1 : 1"<< endl;
      fout << "2 : 2"<< endl;
      fout << "6 : 1 2"<< endl;
      fout << "3 : 3"<< endl;
      fout << "5 : 1 3"<< endl;
      fout << "4 : 2 3"<< endl;
      fout << "#END"<< endl;
    }
  return 1;
}
int main(int argc, char ** argv)
{
  ifstream ale_stream (argv[1]);
  ifstream keep_stream (argv[2]);

  string ale_name=argv[1];
  vector<string> tokens;
  boost::split(tokens,ale_name,boost::is_any_of("."),boost::token_compress_on);
  ale_name=tokens[0];

  string footer=argv[2];
  boost::split(tokens,footer,boost::is_any_of("/"),boost::token_compress_on);
  footer=tokens[tokens.size()-1]; 

  ale_name=ale_name+"_"+footer+".ale";
  
  vector<string> keep_list;
  vector<string> keep_leaves;

  string line;
  while(! keep_stream.eof())
    {      
      getline (keep_stream,line);
      if (line.find("(")!=line.npos)
	{
	  tree_type * T=TreeTemplateTools::parenthesisToTree(line,false,"ID");
	  vector <Node *> leaves=T->getLeaves();
	  for (vector <Node *>::iterator it=leaves.begin();it!=leaves.end();it++)
	    {
	      string name=(*it)->getName();
	      boost::trim(name);
	      vector<string> tokens;
	      boost::split(tokens,name,boost::is_any_of("_"),boost::token_compress_on);
	      string sp=tokens[0];
	      keep_list.push_back(sp);
	    }
	  break;
	}
      else
	{
	  boost::trim(line);
	  vector<string> tokens;
	  boost::split(tokens,line,boost::is_any_of("_"),boost::token_compress_on);
	  string sp=tokens[0];
	  keep_list.push_back(sp);
	}
    }
  
  string field="";
  string constructor_string;
  int N_leaves;
  int observations;
  long int last_leafset_id;
  map <long int,scalar_type> Bip_counts;
  map <long int,map <long int, map <long int,scalar_type > > > Dip_counts;
  map < boost::dynamic_bitset<>,long int>  set_ids;         
  map< long int, boost::dynamic_bitset<> > id_sets;         
  
  map <long int,string> id_leaves;
  map <string,long int> leaf_ids;
  while(! ale_stream.eof())
    {
      getline (ale_stream,line);
      boost::trim(line);
      if (line[0]=='#') field=line;
      else
	{
	  if (field=="#constructor_string")
	    {
	      constructor_string=line;
	    }
	  else if (field=="#observations")
	    {
              boost::trim(line);
              observations=atof(line.c_str());
	    }
	  else if (field=="#Bip_counts")
	    {
	      vector<string> tokens;
	      boost::trim(line);
	      boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
              Bip_counts[atol(tokens[0].c_str())]=atof(tokens[1].c_str());
	    }
	  else if (field=="#Dip_counts")
	    {
	      vector<string> tokens;
              boost::trim(line);
              boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
	      Dip_counts[atol(tokens[0].c_str())][atol(tokens[1].c_str())][atol(tokens[2].c_str())]=atof(tokens[3].c_str());
	    }
	  else if (field=="#last_leafset_id")
	    {
	      boost::trim(line);
	      last_leafset_id=atol(line.c_str());	      
	    }
	  else if (field=="#leaf-id")
	    {
	      vector<string> tokens;
	      boost::trim(line);
	      boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
	      long int id=atol(tokens[1].c_str());
	      string leaf_name=tokens[0];
	      leaf_ids[leaf_name]=id;
	      id_leaves[id]=leaf_name;
	    }
	  else if (field=="#set-id")
	    {
              vector<string> fields;
              boost::trim(line);
              boost::split(fields,line,boost::is_any_of(":"),boost::token_compress_on);
              boost::trim(fields[0]);
              long int set_id=atol(fields[0].c_str());
              vector<string> tokens;
              boost::trim(fields[1]);
              boost::split(tokens,fields[1],boost::is_any_of("\t "),boost::token_compress_on);
              boost::dynamic_bitset<> temp( leaf_ids.size()+1 );
              
              for (vector<string>::iterator it=tokens.begin();it!=tokens.end();it++) { //Setting the proper bits to 1
		temp[static_cast<int>(atoi((*it).c_str()))] = 1; //cout << id_leaves[static_cast<int>(atoi((*it).c_str()))] << " ";
              }

	      //std::cout <<"setid : "<< set_id << " READING: " << temp << std::endl;
              set_ids[temp]=set_id;
              id_sets[set_id]=temp;
	    }
	  
	}
    }
  N_leaves=prune(&constructor_string,keep_list,&keep_leaves);

  if (N_leaves<4) return cout_small(keep_leaves,ale_name);
      

  boost::dynamic_bitset<> keep_set( leaf_ids.size()+1 );
  long int new_last_leafset_id=1;
  map <long int,string> new_id_leaves;
  map <string,long int> new_leaf_ids;
  map <long int,map <long int, map <long int,scalar_type > > > new_Dip_counts;
  map < boost::dynamic_bitset<>,long int>  new_set_ids;         
  map< long int, boost::dynamic_bitset<> > new_id_sets;
  map< long int, long int > new_ids;         

  map < boost::dynamic_bitset<> , boost::dynamic_bitset<> > new_sets;
  map < boost::dynamic_bitset<> ,scalar_type> new_Bip_counts;


  for (vector<string>::iterator it=keep_leaves.begin() ;it!=keep_leaves.end() ; it++)
    {
      string name=(*it);
      boost::dynamic_bitset<> old_set( leaf_ids.size()+1 );
      size_t i = leaf_ids[name];
      keep_set[i] = 1; //! XX !
      old_set[i] = 1; //! XX !

      cout << i << " " << name << " " << leaf_ids[name] <<"  "<< keep_set << " " << old_set << endl;

      new_ids[i]=new_last_leafset_id;
      new_leaf_ids[name] = new_last_leafset_id;
      new_id_leaves[new_last_leafset_id] = name;

      new_set_ids[old_set]=new_last_leafset_id;
      new_id_sets[new_last_leafset_id]=old_set;
      new_sets[old_set] = old_set;

      new_last_leafset_id++;

    }


  for (map < boost::dynamic_bitset<>,long int>::iterator it = set_ids.begin(); it!=set_ids.end(); it++)
    {
      boost::dynamic_bitset<>  old_set = (*it).first;
      boost::dynamic_bitset<>  new_set = old_set & keep_set;
      new_sets[old_set] = new_set;
      //cout <<(*it).second<<" "<< old_set << " -> " << new_set << endl;
      //if (old_set!=keep_set and not (new_set.none())) new_sets[old_set] = new_set;

      if (old_set!=keep_set and new_set!=keep_set and not (new_set.none()) and (new_set_ids.count(new_set)==0) and Bip_counts.count((*it).second))
	new_Bip_counts[new_set] = 0;
    }

  for (map < boost::dynamic_bitset<> ,scalar_type> ::iterator it=new_Bip_counts.begin(); it!=new_Bip_counts.end(); it++)
    {
      new_set_ids[(*it).first]=new_last_leafset_id;
      new_id_sets[new_last_leafset_id]=(*it).first;
      new_last_leafset_id++;
    }
  
  for (map <long int,map <long int, map <long int,scalar_type > > >::iterator it = Dip_counts.begin(); it!=Dip_counts.end(); it++)
    {
      boost::dynamic_bitset<> gamma = id_sets[(*it).first];
      for (map <long int, map <long int,scalar_type > > ::iterator jt = (*it).second.begin(); jt!=(*it).second.end(); jt++)
	{
	  boost::dynamic_bitset<> gamma_p = id_sets[(*jt).first];
	  for ( map <long int,scalar_type >  ::iterator kt = (*jt).second.begin(); kt!=(*jt).second.end(); kt++)
	    {
	      boost::dynamic_bitset<> gamma_pp = id_sets[(*kt).first];

	      if (new_sets[gamma].none() or new_sets[gamma_p].none() or new_sets[gamma_pp].none() or (new_sets[gamma]==keep_set))
		{
		  ;
		}
	      else
		{
		  long int new_g_id   = new_set_ids[new_sets[gamma   ]];
		  long int new_gp_id  = new_set_ids[new_sets[gamma_p ]];
		  long int new_gpp_id = new_set_ids[new_sets[gamma_pp]];
		  //cout << "---------------------------------" << endl;
		  //cout << keep_set << " " << keep_set << " " << keep_set << endl;
		  //cout << gamma << " " << gamma_p << " " << gamma_pp << endl;
		  //cout << new_sets[gamma] << " " << new_sets[gamma_p] << " " << new_sets[gamma_pp] << endl;
		  //cout << "---------------------------------" << endl;
		  new_Bip_counts[new_sets[gamma   ]]+= (*kt).second;
		  if (new_gp_id<new_gpp_id)
		    new_Dip_counts[new_g_id][new_gp_id][new_gpp_id] += (*kt).second;
		  else
		    new_Dip_counts[new_g_id][new_gpp_id][new_gp_id] += (*kt).second;
		}
	    }
	}
    }
  ofstream fout(ale_name.c_str());
  cout << "#writing "<< ale_name << endl;

  fout << "#constructor_string"<< endl;
  fout << constructor_string;

  fout << "#observations"<< endl;
  fout << observations<< endl;

  fout << "#Bip_counts" << endl;
  for (map < boost::dynamic_bitset<> ,scalar_type>::iterator it=new_Bip_counts.begin(); it!=new_Bip_counts.end(); it++)
    {
      fout << new_set_ids[(*it).first] << "\t" << (*it).second << endl;      
    }

  fout << "#Bip_bls"<< endl;
  for (map < boost::dynamic_bitset<>,long int>  ::iterator it=new_set_ids.begin() ;it!=new_set_ids.end() ; it++)         
    {
      long int g_id=(*it).second;
      //fout << (*it).first << "\t";
      fout << g_id << "\t" << 1 << endl;      
    }

  fout << "#Dip_counts"<< endl;
  for (map <long int,map <long int, map <long int,scalar_type > > >::iterator it = new_Dip_counts.begin(); it!=new_Dip_counts.end(); it++)
    {     
      boost::dynamic_bitset<> gamma = new_id_sets[(*it).first];

      for (map <long int, map <long int,scalar_type > > ::iterator jt = (*it).second.begin(); jt!=(*it).second.end(); jt++)
	{
	  boost::dynamic_bitset<> gamma_p = new_id_sets[(*jt).first];

	  for ( map <long int,scalar_type >  ::iterator kt = (*jt).second.begin(); kt!=(*jt).second.end(); kt++)
	    {
	      boost::dynamic_bitset<> gamma_pp = new_id_sets[(*kt).first];
	      long int new_g_id   = new_set_ids[gamma   ];
	      long int new_gp_id  = new_set_ids[gamma_p ];
	      long int new_gpp_id = new_set_ids[gamma_pp];
	      fout << new_g_id << "\t" << new_gp_id << "\t" << new_gpp_id << "\t" << (*kt).second << endl;  	      
	    }
	}
    }
  
  fout << "#last_leafset_id"<< endl;
  fout << new_last_leafset_id-1 << endl;

  fout << "#leaf-id"<< endl;
  for (vector<string>::iterator it=keep_leaves.begin() ;it!=keep_leaves.end() ; it++)
    {
      fout << (*it) << "\t" <<  new_leaf_ids[(*it)] <<endl; 
    }

  fout << "#set-id"<< endl;
  for (map < boost::dynamic_bitset<>,long int>  ::iterator it=new_set_ids.begin() ;it!=new_set_ids.end() ; it++)         
    {
      boost::dynamic_bitset<> gamma=(*it).first;
      long int g_id=(*it).second;
      //fout << gamma << "\t:";
      fout << g_id << "\t:";
      
      size_t i = gamma.find_first();
      while(i != boost::dynamic_bitset<>::npos)
	{
	  fout << "\t" << new_ids[i];
	  i = gamma.find_next(i);
	}
      fout << endl;
    }
  fout << "#END"<< endl;
  return 1;


  
  
}
