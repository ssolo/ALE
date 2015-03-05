#include "mpi_tree.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
using namespace boost::mpi;

void mpi_tree::load_distributed_ales(string fname)
{
  client_fnames.clear();
  vector<vector<string> > scatter_fnames;//del-loc
  if (rank==server)
    {
      ifstream file_stream (fname.c_str());
      int tree_i=0;
      set <string> verify;
      if (file_stream.is_open())  //  ########## read trees ############
	{
	  while (! file_stream.eof())
	    {
	      string line;
	      getline (file_stream,line);
	      vector <string> tokens;
	      boost::trim(line);	    
	      boost::split(tokens,line,boost::is_any_of("\t "),boost::token_compress_on);
	      int client_i=atoi(tokens[1].c_str());
	      string ale_file=tokens[0];
	      if (ale_file.size()>1)
		{
		  if (not (client_i<scatter_fnames.size())) {vector<string> tmp; scatter_fnames.push_back(tmp);}
		  scatter_fnames[client_i].push_back(ale_file);
		  verify.insert(ale_file);
		}
	    }
	}
      cout << "# Scattering: " << verify.size() << " ale files.."<<endl;
      N_ales=verify.size() ;
      verify.clear();
    }
  scatter(world,scatter_fnames,client_fnames,server);

  if (rank==server) cout << "#..loading.." << endl;
  
  for ( vector<string>::iterator it=client_fnames.begin();it!=client_fnames.end();it++)
    {
      //cout << rank << " has " << (*it) << endl;
      approx_posterior * ale;//del-loc
      ale = load_ALE_from_file((*it));      
      ale_pointers.push_back(ale);      
    }

  //del-locs
  for ( vector<vector<string> >::iterator jt=scatter_fnames.begin();jt!=scatter_fnames.end();jt++) 
    (*jt).clear();

  scatter_fnames.clear();

  if (rank==server) cout << "# done." <<endl;
  scalar_type tmp=N_ales;
  broadcast(world,tmp,server);
  N_ales=tmp;


}
void mpi_tree::distribute_ales(vector<string> fnames,bool list_of_trees)
{

  client_fnames.clear();
  vector<vector<string> > scatter_fnames;//del-loc
  
  if (rank==server)
    {
      //cout << "#rank:" <<rank <<" of " << size <<" started distribute."<<endl;
      
      set <string> verify;
      for (vector<string>::iterator it=fnames.begin();it!=fnames.end();it++)
	verify.insert((*it));
      cout << "# Distributing: " << verify.size() << " ale files.."<<endl;
      for (int i=0;i<size;i++)
	{
	  vector<string> tmp;
	  scatter_fnames.push_back(tmp);
	}
      map <string,int> fname_counts;//del-loc
      map <int,vector<string> > count_fnames;//del-loc
      vector<string> sorted_fnames;//del-loc
      
      for (vector<string>::iterator it=fnames.begin();it!=fnames.end();it++)
	{
	  string fname=(*it);
	  approx_posterior * ale;
	  if (not list_of_trees or fname.find(".ale")!=fname.npos) 
	    ale = load_ALE_from_file(fname);      
	  else 
	    ale = observe_ALE_from_string(fname);      
	  int gid_count=0;
	  for (int i=0;i<(int)ale->Dip_counts.size();i++)
	    gid_count+=ale->Dip_counts[i].size();
	  gid_count+=ale->Bip_counts.size();
	  if (gid_count<1e5)//we could be carful to fiter big ales ... probably should automate this top 1%?
	    {
	      fname_counts[fname]=gid_count;
	      count_fnames[gid_count].push_back(fname);
	    }
	  delete ale;
	}
      for (map <int,vector<string> >::iterator jt=count_fnames.begin();jt!=count_fnames.end();jt++) 
	for (vector<string> ::iterator kt=(*jt).second.begin();kt!=(*jt).second.end();kt++) 	  
	  sorted_fnames.push_back((*kt));
      for (int i=0;i<(int)sorted_fnames.size();i++)
	scatter_fnames[i%size].push_back(sorted_fnames[i]);
      

      map <scalar_type,int> gidsum_ranks;//del-loc
      //we try to exhange fnames to optimze gid distribution
      while (1)
	{
	  int max_rank=-1;
	  int min_rank=-1;
	  scalar_type max_sum=0; 
	  scalar_type min_sum=6e23; 

	  for (int i=0;i<size;i++)
	    {
	      scalar_type gidsum=0;
	      for (int j=0;j<(int)scatter_fnames[i].size();j++)
		gidsum+=fname_counts[scatter_fnames[i][j]];
	      while (gidsum_ranks.count(gidsum)!=0)
		gidsum+=0.1;
	      gidsum_ranks[gidsum]=i; 
	      //cout << i << " has " << gidsum << " " <<scatter_fnames[i].size()  << endl;
	      if (max_sum<gidsum)
		{
		  max_sum=gidsum;
		  max_rank=i;
		}
	      if (min_sum>gidsum)
		{
		  min_sum=gidsum;
		  min_rank=i;
		}
	    }
	  //cout << endl;
	  max_sum=0;
	  for (int j=0;j<(int)scatter_fnames[max_rank].size();j++)
	    max_sum+=fname_counts[scatter_fnames[max_rank][j]];
	  min_sum=0;
	  for (int j=0;j<(int)scatter_fnames[min_rank].size();j++)
	    min_sum+=fname_counts[scatter_fnames[min_rank][j]];
	  //cout << max_rank << " and " << min_rank << endl;
	  bool changed=false;
	  for (int j=0;j<(int)scatter_fnames[max_rank].size();j++)
	    for (int k=0;k<(int)scatter_fnames[min_rank].size();k++)
	      {
		scalar_type max_rank_count=fname_counts[scatter_fnames[max_rank][j]];
		scalar_type min_rank_count=fname_counts[scatter_fnames[min_rank][k]];
		if (abs(min_sum-max_sum) > abs( (min_sum-min_rank_count+max_rank_count) - (max_sum+min_rank_count-max_rank_count ) ))
		  {
		    string jname=scatter_fnames[max_rank][j];
		    string kname=scatter_fnames[min_rank][k];
		    scatter_fnames[max_rank].insert(scatter_fnames[max_rank].begin()+j,kname);
		    scatter_fnames[min_rank].insert(scatter_fnames[min_rank].begin()+k,jname);
		    scatter_fnames[max_rank].erase(scatter_fnames[max_rank].begin()+j+1);
		    scatter_fnames[min_rank].erase(scatter_fnames[min_rank].begin()+k+1);

		    min_sum=min_sum-min_rank_count+max_rank_count;
		    max_sum=max_sum+min_rank_count-max_rank_count;
		    changed=true;
		  }
	      }
	  for (int j=0;j<(int)scatter_fnames[max_rank].size();j++)
	    {
	      scalar_type max_rank_count=fname_counts[scatter_fnames[max_rank][j]];
	      if (abs(min_sum-max_sum) > abs( (min_sum+max_rank_count) - (max_sum-max_rank_count ) ))
		{
		  scatter_fnames[min_rank].push_back(scatter_fnames[max_rank][j]);
		  scatter_fnames[max_rank].erase(scatter_fnames[max_rank].begin()+j);		  
		  changed=true;		  
		  break;
		}
	    }
	  gidsum_ranks.clear();
	  if (not changed) break;	    
	}

      verify.clear();
      for (int i=0;i<size;i++)
	for (int j=0;j<(int)scatter_fnames[i].size();j++)
	  verify.insert(scatter_fnames[i][j]);
      //del-locs

      fname_counts.clear();
      for (map <int,vector<string> >::iterator jt=count_fnames.begin();jt!=count_fnames.end();jt++) 
	(*jt).second.clear();
      count_fnames.clear();
      sorted_fnames.clear();
      gidsum_ranks.clear();

      cout << "# Scattering: " << verify.size() << " ale files.."<<endl;
      N_ales=verify.size() ;
      verify.clear();
    }

  scatter(world,scatter_fnames,client_fnames,server);

  if (rank==server) cout << "#..loading.." << endl;
  
  for ( vector<string>::iterator it=client_fnames.begin();it!=client_fnames.end();it++)
    {
      approx_posterior * ale;//del-loc
      if (not list_of_trees or (*it).find(".ale")!=(*it).npos) 
	ale = load_ALE_from_file((*it));      
      else 
	ale = observe_ALE_from_string((*it));      

      if (scalar_parameter["use_mpp_trees"]==0)
	{
	  ale_pointers.push_back(ale);      
	}
      else
	{
	  cout<< "Using mpp trees..!!" << endl;
	  vector <string >trees;
	  trees.push_back(ale->mpp_tree().first);
	  ale_pointers.push_back(observe_ALE_from_strings(trees));      
	  delete ale;
	}
    }

  //del-locs
  for ( vector<vector<string> >::iterator jt=scatter_fnames.begin();jt!=scatter_fnames.end();jt++) 
    (*jt).clear();


  scatter_fnames.clear();

  if (rank==server) cout << "# done." <<endl;
  scalar_type tmp=N_ales;
  broadcast(world,tmp,server);
  N_ales=tmp;
}
void mpi_tree::clear_counts()
{
  for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
    {
      string count_name=(*it).first;
      for (int branch=0;branch<model->last_branch;branch++)	
	model->branch_counts[count_name][branch]=0;
    }
  for (map<string,scalar_type> ::iterator it=model->MLRec_events.begin();it!=model->MLRec_events.end();it++)
    {
      model->MLRec_events[(*it).first]=0;
    }

  model->Ttokens.clear();

  int done=1;
  broadcast(world,done,server);
}
void mpi_tree::gather_T_to_from()
{
  if (rank==server)
    {      
      gathered_T_to_from.clear();
      for (int e=0;e<model->last_branch;e++)
	{
	  vector <scalar_type> tmp;
	  gathered_T_to_from.push_back(tmp);
	  for (int f=0;f<model->last_branch;f++)
	    gathered_T_to_from[e].push_back(0);
	}
    }
  vector< vector < vector<scalar_type> > >  gather_vector;
  gather(world,model->T_to_from,gather_vector,server);
  if (rank==server)
    {
      scalar_type Tsum=0;
      for (vector< vector < vector<scalar_type> > > ::iterator it=gather_vector.begin();it!=gather_vector.end();it++)
	for (int e=0;e<model->last_branch;e++)
	  for (int f=0;f<model->last_branch;f++)
	    {
	      gathered_T_to_from[e][f]+=(*it)[e][f];
	      Tsum+=(*it)[e][f];
	    }
      cout << "TOTAL Ts = "<< Tsum << endl;
      //map <scalar_type, vector< int > >sort_e;
      //map <scalar_type, vector< int > >sort_f;
      sort_e.clear();
      sort_f.clear();
      for (int e=0;e<model->last_branch;e++)
	for (int f=0;f<model->last_branch;f++)
	  {
	    sort_e[-gathered_T_to_from[e][f]].push_back(e);
	    sort_f[-gathered_T_to_from[e][f]].push_back(f);	      
	  }
      for (map <scalar_type, vector< int > >::iterator it=sort_e.begin();it!=sort_e.end();it++)
	{
	  scalar_type Ts=(*it).first;
	  for (int i=0;i<(*it).second.size();i++)
	    {
	      int e=sort_e[Ts][i];
	      int f=sort_f[Ts][i];
	      if (Ts< -0.01*Tsum*0.1 and false)
		{
		  cout << Ts;	      
		  if (e<model->last_leaf)
		    cout << " " << model->node_name[model->id_nodes[e]];
		  else
		    cout << " " << e;
		  if (f<model->last_leaf)
		    cout << " " << model->node_name[model->id_nodes[f]];	      
		  else
		    cout << " " << f;
		  cout << endl;
		}
	    }
	}
    }
}
void mpi_tree::gather_counts()
{
  map< string ,vector <vector <scalar_type> > > gathered_branch_counts;//del-loc
  for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
    {
      string count_name=(*it).first;
      if (rank==server)
	{
	  vector <vector <scalar_type> > tmp;
	  gathered_branch_counts[count_name]=tmp;
	}
      gather(world,model->branch_counts[count_name],gathered_branch_counts[count_name],server);
      if (rank==server)
	{
	  for (int branch=0;branch<model->last_branch;branch++)	
	    {
	      for (int i=1;i<size;i++)
		model->branch_counts[count_name][branch]+=gathered_branch_counts[count_name][i][branch];
	    }    
	  model->show_counts(count_name);
	}  
    }  
  //del-locs

  for ( map< string ,vector <vector <scalar_type> > >::iterator it= gathered_branch_counts.begin();it!= gathered_branch_counts.end();it++)
    {
      for ( vector <vector <scalar_type> >::iterator jt= (*it).second.begin();jt!= (*it).second.end();jt++)
	(*jt).clear();
      (*it).second.clear();
    }
  gathered_branch_counts.clear();

  //Ttokens
  vector<vector<string> > gather_Ttokens;  //del-loc
  gather(world,model->Ttokens,gather_Ttokens,server);
  if (rank==server) for (int i=0;i<size;i++) for (int j=0;j<(int)gather_Ttokens[i].size();j++) Ttokens.push_back(gather_Ttokens[i][j]);

  //del-locs
 
  for (vector<vector<string> >::iterator it=gather_Ttokens.begin();it!=gather_Ttokens.end();it++)
    (*it).clear();
  gather_Ttokens.clear(); 



  int done=1;
  broadcast(world,done,server);
}


void mpi_tree::print_branch_counts()
{
  if (rank==server)
    {
      cout<<"#\t";
      for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
	{
	  string count_name=(*it).first;
	  cout << count_name << "\t";
	}
      cout << endl;
      for (int branch=0;branch<model->last_branch;branch++)
	{
	  cout << branch << "\t";
	  for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
	    {
	      string count_name=(*it).first;
	      cout << model->branch_counts[count_name][branch] << "\t";
	    }
	  cout << endl;
	}
    }
}
string mpi_tree::branch_counts_string()
{
  string return_string;
  if (rank==server)
    {
      stringstream out; 

      scalar_type total_D=0,total_T=0,total_L=0,total_S=0;
      for (int branch=0;branch<model->last_branch;branch++)
	{
	  total_D+=model->branch_counts["Ds"][branch];
	  total_T+=model->branch_counts["Ts"][branch];
	  total_L+=model->branch_counts["Ls"][branch];
	  total_S+=model->branch_counts["copies"][branch];	  
	}
      out << " " <<total_D;
      out << " " <<total_T;
      out << " " <<total_L;
      out << " " <<total_S;
      return_string=out.str();
    }
  broadcast(world,return_string,server);

  return return_string;
}
void mpi_tree::show_branch_counts()
{
  if (rank==server)
    {
      for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
	{
	  string count_name=(*it).first;
	  model->show_counts(count_name);
	}
      scalar_type total_D=0,total_T=0,total_L=0,total_S=0;
      for (int branch=0;branch<model->last_branch;branch++)
	{
	  total_D+=model->branch_counts["Ds"][branch];
	  total_T+=model->branch_counts["Ts"][branch];
	  total_L+=model->branch_counts["Ls"][branch];
	  total_S+=model->branch_counts["copies"][branch];

	}
      cout << "#total D: " << total_D;
      cout << " T: " << total_T;
      cout << " L: " <<total_L;
      cout << " S: " <<total_S;

      cout << endl;
      cout << "#avg. D: " << total_D/N_ales;
      cout << " T: " << total_T/N_ales;
      cout << " L: " << total_L/N_ales;
      cout << " S: " <<total_S/N_ales/((model->last_branch+1));
      
      cout << endl;      
    }
  int done=1;
  broadcast(world,done,server);
}
scalar_type mpi_tree::calculate_MLRecs(bool estimate,bool branchwise)
{
  clear_counts();
  MLRec_res.clear();
  scalar_type ll=0;
  vector<scalar_type> gather_ll;
  //show_rates();
  boost::timer * t = new boost::timer();
  for (int i=0;i<(int)ale_pointers.size();i++)
    {
      //cout << rank <<" at " <<round(i/(float)ale_pointers.size()*100.)<<" %, strats "<< client_fnames[i] << endl;
      pair <string,scalar_type > res=model->p_MLRec(ale_pointers[i]);
      if (model->signal==-11)
	{
	  model->signal=0;
	  cout << "ERERERER "<<rank << " " << i << endl;
	  ale_pointers[i]->save_state("error_ale");
	}
      cout <<"#ML " << client_fnames[i] << " " << res.first << " " << log(res.second) << endl;
      ll+=log(res.second);
      MLRec_res.push_back(res);      
    }

 cout << rank << " "<<t->elapsed() <<endl;

  
  
  gather(world,ll,gather_ll,server);
  scalar_type sum_ll=0;
  if (rank==server) for (int i=0;i<size;i++) sum_ll+=gather_ll[i];
  broadcast(world,sum_ll,server);

  gather_counts();

  return sum_ll;
}

scalar_type mpi_tree::calculate_p()
{  
  scalar_type ll=0;
  model->calculate_EGb();
  vector<scalar_type> gather_ll;

  boost::timer * t = new boost::timer();
  for (int i=0;i<(int)ale_pointers.size();i++)
    {
      //cout << rank <<" at " <<round(i/(float)ale_pointers.size()*100.)<<" %, strats "<< client_fnames[i] << endl;
      scalar_type tmpp=model->p(ale_pointers[i]);
      if (tmpp==0) cout << client_fnames[i] << " is 0 !!"<<endl;
      cout <<"#LL " << client_fnames[i] << " " << log(tmpp) << endl;

      ll +=log(tmpp);
    }
  //cout << rank << " "<<t->elapsed() <<endl;
  gather(world,ll,gather_ll,server);

  scalar_type sum_ll=0;
  if (rank==server) for (int i=0;i<size;i++) sum_ll+=gather_ll[i];
  broadcast(world,sum_ll,server);
  if (rank==server) cout << model->scalar_parameter["delta_avg"] <<  " " << model->vector_parameter["N"][0]*model->scalar_parameter["tau_avg"] << " " << model->scalar_parameter["lambda_avg"] << " " << model->vector_parameter["Delta_bar"][0] <<" " << sum_ll<<endl;

  return sum_ll;
}
void mpi_tree::estimate_rates()
{
  scalar_type delta,tau,lambda;
  map< string ,vector <vector <scalar_type> > > gathered_branch_counts;//del-loc
  for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
    {
      string count_name=(*it).first;
      if (rank==server)
	{
	  vector <vector <scalar_type> > tmp;
	  gathered_branch_counts[count_name]=tmp;
	}
      gather(world,model->branch_counts[count_name],gathered_branch_counts[count_name],server);
      if (rank==server)
	{
	  for (int branch=0;branch<model->last_branch;branch++)	
	    {
	      for (int i=1;i<size;i++)
		model->branch_counts[count_name][branch]+=gathered_branch_counts[count_name][i][branch];
	    }    
	  //model->show_counts(count_name);
	}  
    }  
  if (rank==server)
    {
      
      scalar_type Csum=0;
      for (int e=0;e<model->last_branch;e++)	
	{
	  Csum+=model->branch_counts["count"][e];
	}
      scalar_type P_D_avg=0;
      scalar_type P_T_avg=0;
      scalar_type P_L_avg=0;
      scalar_type w_sum=0;
  
      for (int e=0;e<model->last_branch;e++)	
	{
	  scalar_type N_S=model->branch_counts["count"][e];
	  scalar_type Ee=model->branch_counts["Ls"][e]/N_S;
	  scalar_type Ge=model->branch_counts["singleton"][e]/N_S;
	  //cout << Ee << " " << Ge << endl;
	  scalar_type P_D=max(1e-6,-1*((1 - Ee - Ge)/((-1 + Ee)*(-2*Ee + Ge + Ee*Ge))));
	  scalar_type P_L=max(1e-6,-1*((Ee - 3*Ee*Ee + Ee*Ee*Ge)/((-1 + Ee)*(-2*Ee + Ge + Ee*Ge))));	    
	  scalar_type P_T=max(1e-6,model->branch_counts["Ts"][e]/Csum*(float)model->last_branch);
	  //cout << e<<" "<< P_D << " " << P_L << " " << P_T << " " << " "<< model->branch_counts["count"][e] << endl;
	  scalar_type w=1;//model->branch_counts["count"][e];
	  P_D_avg+=( P_D ) * w;
	  P_T_avg+=( P_T ) * w;
	  P_L_avg+=( P_L ) * w;
	  w_sum+=w;
	}
      delta=P_D_avg/w_sum;
      tau=P_T_avg/w_sum;
      lambda=P_L_avg/w_sum;
      cout << " rate estimates " << delta << " " << tau << " " << lambda << endl;
    }
  broadcast(world,delta,server);
  broadcast(world,tau,server);
  broadcast(world,lambda,server);  
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau",tau);
  model->set_model_parameter("lambda",lambda);
}

void mpi_tree::estimate_rates_bw()
{
  vector<scalar_type> delta,tau,lambda;
  map< string ,vector <vector <scalar_type> > > gathered_branch_counts;//del-loc
  for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
    {
      string count_name=(*it).first;
      if (rank==server)
	{
	  vector <vector <scalar_type> > tmp;
	  gathered_branch_counts[count_name]=tmp;
	}
      gather(world,model->branch_counts[count_name],gathered_branch_counts[count_name],server);
      if (rank==server)
	{
	  for (int branch=0;branch<model->last_branch;branch++)	
	    {
	      for (int i=1;i<size;i++)
		model->branch_counts[count_name][branch]+=gathered_branch_counts[count_name][i][branch];
	    }    
	  //model->show_counts(count_name);
	}  
    }  
  if (rank==server)
    {
      
      scalar_type Csum=0;
      for (int e=0;e<model->last_branch;e++)	
	{
	  Csum+=model->branch_counts["count"][e];
	}
      scalar_type P_D_avg=0;
      scalar_type P_T_avg=0;
      scalar_type P_L_avg=0;
      scalar_type w_sum=0;
  
      for (int e=0;e<model->last_branch;e++)	
	{
	  scalar_type N_S=model->branch_counts["count"][e];
	  scalar_type Ee=model->branch_counts["Ls"][e]/N_S;
	  scalar_type Ge=model->branch_counts["singleton"][e]/N_S;
	  //cout << Ee << " " << Ge << endl;
	  scalar_type P_D=max(1e-6,-1*((1 - Ee - Ge)/((-1 + Ee)*(-2*Ee + Ge + Ee*Ge))));
	  scalar_type P_L=max(1e-6,-1*((Ee - 3*Ee*Ee + Ee*Ee*Ge)/((-1 + Ee)*(-2*Ee + Ge + Ee*Ge))));	    
	  scalar_type P_T=max(1e-6,model->branch_counts["Ts"][e]/Csum*(float)model->last_branch);
	  delta.push_back( P_D );
	  tau.push_back( P_T);
	  lambda.push_back( P_L );
	  cout << " rate estimates " << e <<" " << delta[e] << " " << tau[e] << " " << lambda[e] << endl;
	}
    }
  broadcast(world,delta,server);
  broadcast(world,tau,server);
  broadcast(world,lambda,server);  
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau",tau);
  model->set_model_parameter("lambda",lambda);
}

scalar_type mpi_tree::calculate_pun(int n, bool bw)
{
  scalar_type ll=calculate_pun();
  
  for (int i=0;i<n;i++)
    {
      if (bw) estimate_rates_bw(); else estimate_rates();
      ll=calculate_pun();
      if (rank==server) cout << ll << " " << bw << endl;
    }
  return ll;
}
scalar_type mpi_tree::calculate_pun()
{  
  scalar_type ll=0;
  model->calculate_undatedEs();
  vector<scalar_type> gather_ll;
  model->MLRec_events.clear();
  for (int e=0;e<model->last_branch;e++)	
    {
      model->branch_counts["Os"][e]=0;
      model->branch_counts["Ds"][e]=0;
      model->branch_counts["Ts"][e]=0;
      model->branch_counts["Tfroms"][e]=0;     
      model->branch_counts["Ls"][e]=0;
      model->branch_counts["count"][e]=0;
      model->branch_counts["copies"][e]=0;
      model->branch_counts["singleton"][e]=0;
    }
  for (int e=0;e<model->last_branch;e++)
      for (int f=0;f<model->last_branch;f++)
	model->T_to_from[e][f]=0;

  boost::timer * t = new boost::timer();
  for (int i=0;i<(int)ale_pointers.size();i++)
    {
      //cout << rank <<" at " <<round(i/(float)ale_pointers.size()*100.)<<" %, strats "<< client_fnames[i] << endl;
      scalar_type tmpp=model->pun(ale_pointers[i]);
      for (int i=0;i<100;i++) model->sample_undated();
      if (tmpp==0) cout << client_fnames[i] << " is 0 !!"<<endl;
      //cout <<"#LL " << client_fnames[i] << " " << log(tmpp) << endl;
      
      ll +=log(tmpp);
    }
  //cout << rank << " "<<t->elapsed() <<endl;
  gather(world,ll,gather_ll,server);

  scalar_type sum_ll=0;
  if (rank==server) for (int i=0;i<size;i++) sum_ll+=gather_ll[i];
  broadcast(world,sum_ll,server);
  //if (rank==server) cout <<" " << sum_ll<<endl;

  return sum_ll;
}

scalar_type mpi_tree::calculate_punt(string S)
{  
  scalar_type ll=0;
  model->calculate_undatedEs();
  vector<scalar_type> gather_ll;
  model->MLRec_events.clear();
  for (int e=0;e<model->last_branch;e++)	
    {
      model->branch_counts["Os"][e]=0;
      model->branch_counts["Ds"][e]=0;
      model->branch_counts["Ts"][e]=0;
      model->branch_counts["Tfroms"][e]=0;     
      model->branch_counts["Ls"][e]=0;
      model->branch_counts["count"][e]=0;
      model->branch_counts["copies"][e]=0;
    }
  for (int e=0;e<model->last_branch;e++)
      for (int f=0;f<model->last_branch;f++)
	model->T_to_from[e][f]=0;

  boost::timer * t = new boost::timer();
  for (int i=0;i<(int)ale_pointers.size();i++)
    {
      //cout << rank <<" at " <<round(i/(float)ale_pointers.size()*100.)<<" %, strats "<< client_fnames[i] << endl;
      exODT_model * alt_model=new exODT_model();
      alt_model->construct_undated(S);//del-loc
      alt_model->set_model_parameter("delta",scalar_parameter["inital_delta"]);
      alt_model->set_model_parameter("tau",scalar_parameter["inital_tau"]);
      alt_model->set_model_parameter("lambda",scalar_parameter["inital_lambda"]);
      alt_model->calculate_undatedEs();
      scalar_type tmpp=alt_model->pun(ale_pointers[i]);
      for (int i=0;i<10;i++) model->sample_undated();
      if (tmpp==0) cout << client_fnames[i] << " is 0 !!"<<endl;
      //cout <<"#LL " << client_fnames[i] << " " << log(tmpp) << endl;
      
      ll +=log(tmpp);
    }
  //cout << rank << " "<<t->elapsed() <<endl;
  gather(world,ll,gather_ll,server);

  scalar_type sum_ll=0;
  if (rank==server) for (int i=0;i<size;i++) sum_ll+=gather_ll[i];
  broadcast(world,sum_ll,server);
  //if (rank==server) cout <<" " << sum_ll<<endl;

  return sum_ll;
}
