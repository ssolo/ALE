#include "exODT.h"
#include <Bpp/Numeric/Random/RandomTools.h>
using namespace std;
using namespace bpp;

unsigned int good_seed()
{
    unsigned int random_seed, random_seed_a, random_seed_b; 
    std::ifstream file ("/dev/random", std::ios::binary);
    if (file.is_open())
    {
        char * memblock;
        int size = sizeof(int);
        memblock = new char [size];
        file.read (memblock, size);
        file.close();
        random_seed_a = long(memblock);
        delete[] memblock;
    }// end if
    else
    {
        random_seed_a = 0;
    }
    random_seed_b = std::time(0);
    random_seed = random_seed_a xor random_seed_b;
    return random_seed;
} // end good_seed()


int main(int argc, char ** argv)
{

  //simulate
  long int S_seed=good_seed();
  if (argc>1) S_seed=atol(argv[1]);
  long int G_seed=good_seed();
  if (argc>2) G_seed=atol(argv[2]);
  
  int N=1000;
  int n=36;
  int G_n=10;

  scalar_type omega=N*0.08;
  //O rate overall 
  scalar_type delta=0.2;
  //D rate per gene
  scalar_type tau=0.8;
  //T rate per gene
  scalar_type lambda=1.4;
  //L rate per gene

  scalar_type init_t=2;
  scalar_type sigma=N;

  cout << "# Species seed is : " << S_seed << endl;
  cout << "# Genes seed is : "  << G_seed << endl;

  RandomTools::setSeed(S_seed);
  
  vector <long int> population;
  long int next_index=0;
  for (int i=0;i<N;i++)
    {
      population.push_back(i);
      next_index=i;
    }
  next_index++;

  //we record history 
  vector< vector< long int > > families;
  map <long long,scalar_type > event_times;

  vector < int > births;
  vector < int > deaths;


  scalar_type t=init_t;
  long long species_event=0;
  
  int Ds=0;
  int Ts=0;
  int Ls=0;
  int Os=0;
  int Ss=0;

  while (1)
    {
      scalar_type t_next=RandomTools::randExponential(1./(sigma*N));
      t-=t_next;
      if (t<0) break;

      int death=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(N);    
      deaths.push_back(death);

      int birth=death;
      while (birth==death) birth=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(N);    
      births.push_back(birth);      

      vector <long int> family;
      long int mother=population[birth];
      long int daugther=next_index;
      next_index++;
      long int son=next_index;
      next_index++;
      family.push_back(mother);
      family.push_back(daugther);
      family.push_back(son);

      population[birth]=daugther;
      population[death]=son;

      families.push_back(family);
      event_times[species_event]=t;

      species_event++;

    }
  long long number_of_species_events=species_event;
  //sample

  vector <int> population_indicies;
  for (int i=0;i<N;i++) population_indicies.push_back(i);

  vector <int> sampled_population_indicies;  
  for (int i=0;i<n;i++) sampled_population_indicies.push_back(-1);
    
  RandomTools::getSample(population_indicies,sampled_population_indicies);

  vector<long int> sampled_population;
  for (int i=0;i<n;i++) sampled_population.push_back(population[ sampled_population_indicies[i] ]);  

  //traceback

  map<long int,int> sampled_population_counts;
  map<long int,string> strings;
  map<long int,scalar_type> age;
  for (int i=0;i<n;i++) 
    {
      long int extant_species = sampled_population[i];
      stringstream extant_species_name;
      extant_species_name << i;//extant_species;
      strings[extant_species]=extant_species_name.str();
      sampled_population_counts[extant_species]=1;  
      age[extant_species]=0;
    }
  long int lca;
  int rank=0;
  for(vector<vector<long int> >::reverse_iterator event=families.rbegin();event!=families.rend();event++)
    {
      scalar_type t_event=event_times[species_event-1];
      long int mother=(*event)[0];
      long int daugther=(*event)[1];
      long int son=(*event)[2];	      
      if (sampled_population_counts[daugther]==1 and sampled_population_counts[son]==1)
	{
	  rank++;
	  sampled_population_counts[daugther]=0; 
	  sampled_population_counts[son]=0; 
	  sampled_population_counts[mother]=1; 
	  stringstream sons_bl;
	  stringstream daugthers_bl;
	  stringstream rank_bs;
	  rank_bs << rank;
	  sons_bl << t_event - age[son];
	  daugthers_bl << t_event - age[daugther];
	  strings[mother]="("+strings[daugther]+":"+daugthers_bl.str()+","+strings[son]+":"+sons_bl.str()+")"+rank_bs.str();
	  age[mother]=t_event;
	  lca=mother;
	}
      else if (sampled_population_counts[daugther]==1)
	{
	  sampled_population_counts[daugther]=0; 
	  sampled_population_counts[mother]=1; 
	  strings[mother]=strings[daugther];
	  age[mother]=age[daugther];

	}
      else if (sampled_population_counts[son]==1)
	{
	  sampled_population_counts[son]=0; 
	  sampled_population_counts[mother]=1; 
	  strings[mother]=strings[son];
	  age[mother]=age[son];	  
	}	
      species_event--;
    }
  stringstream root_bl;
  root_bl<<age[lca];//init_t-age[lca]; //more general printing anyway needed for genes..

  stringstream fname;
  fname << "S_" << S_seed << ".tree";
  ofstream s_out( fname.str().c_str() );

  s_out  << strings[lca] <<":"<<root_bl.str()<<";"<<endl;
  cout  << strings[lca] <<":"<<root_bl.str()<<";"<<endl;
  
  // rates in units corresponding to unit height S
  //delta  /= age[lca];
  //tau    /= age[lca];
  //lambda /= age[lca];
  
  // a random tree on the same labels
  stringstream rfname;
  rfname << "R_" << S_seed << ".tree";
  ofstream r_out( rfname.str().c_str() );
  vector<string> random_tree_population; 
  map<string,scalar_type> random_tree_ages;
  for (int i=0;i<n;i++) 
  {
    stringstream tmp;
    tmp<<i;
    random_tree_population.push_back(tmp.str());
    random_tree_ages[tmp.str()]=0;
  }
  t=0;
  while(random_tree_population.size()>1)
    {
      int Nr=random_tree_population.size();
      scalar_type t_next=RandomTools::randExponential(1./(2*Nr));
      t+=t_next;
      int i=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(Nr);    
      int j=i;
      while (i==j) j=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(Nr);
      stringstream tmp;
      tmp<<"("<<random_tree_population[i]<<":"<<t-random_tree_ages[random_tree_population[i]]<<","<<random_tree_population[j]<<":"<<t-random_tree_ages[random_tree_population[j]]<<")";       
      random_tree_ages[tmp.str()]=t;
       if (j>i)
	 {
	   random_tree_population.erase(random_tree_population.begin()+j);
	   random_tree_population.erase(random_tree_population.begin()+i);
	 }
       else
	 {
	   random_tree_population.erase(random_tree_population.begin()+i);
	   random_tree_population.erase(random_tree_population.begin()+j);
	 }
       random_tree_population.push_back(tmp.str());	 
    }
  r_out << random_tree_population[0]<<";" << endl;

  // genes
  RandomTools::setSeed(G_seed);

  //we seed genes .. this part is trvilally parallelizable..
  long int gene_count=0;
  map <int,vector<long int> > population_of_genes;
  for (int i=0;i<N;i++) 
    {
      vector <long int> genes_in_species_i;
      //we could have more genes per species .. 
      for (int j=0;j<G_n;j++)
	{	
	  genes_in_species_i.push_back(gene_count);
	  gene_count++;
	}
      population_of_genes[i]=genes_in_species_i;
    };
  
  //..we replay species history 
  species_event=0;
  //and record gene stories..
  long int next_gene=gene_count;
  vector<vector< vector<long int> > > gene_families;
  map <long long, scalar_type > gene_event_times;
  map <long long, string > gene_event_types;

  t=init_t;
  species_event=0;
  cout << "#gene stories .." << endl;
  long long gene_event=0;
  boost::progress_display show_progress( number_of_species_events );
  while(species_event<number_of_species_events)
    {
     
      scalar_type rate_sum=gene_count*(delta+tau+lambda)+omega;
      scalar_type t_next=RandomTools::randExponential(1./(rate_sum));

      scalar_type t_species_event=event_times[species_event];

      //gene event happens
      if (t-t_next > t_species_event)
	{
	  t-=t_next;

	  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(rate_sum);
	  if (r<gene_count*delta+gene_count*tau+gene_count*lambda)
	    {
	      int gene=-1;
	      int species=-1;

	      long int gene_r=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(gene_count);
	      long int gene_re_count=0;
	      gene_re_count=0;
	      for (int i=0;i<N;i++)
		{
		  if (gene_r<gene_re_count+(long int)population_of_genes[i].size())
		    {
		      gene= (gene_re_count+(long int)population_of_genes[i].size()-gene_r) - 1;
		      species=i;
		      break;
		    }		  
		  gene_re_count+=population_of_genes[i].size();
		}
	      
	      if (r<gene_count*delta)
		//D
		{
		  // cout << "D"<<endl;
		  Ds+=1;
		  vector <long int> family;
		  long int father=population_of_genes[species][gene];
		  family.push_back(father);
		  long int daugther=next_gene;
		  next_gene++;
		  family.push_back(daugther);
		  long int son=next_gene;
		  next_gene++;
		  family.push_back(son);

		  vector < vector <long int> > fam_vec;
		  fam_vec.push_back(family);
		  gene_families.push_back(fam_vec);
		  gene_event++;

		  gene_event_times[father]=t;
		  gene_event_types[father]="D";
      
		  population_of_genes[species][gene]=daugther;
		  population_of_genes[species].push_back(son);
		  gene_count++;
		  // cout << "D."<<endl;

		}
	      else if (r<gene_count*delta+gene_count*tau)
		//T
		{		  
		  //cout << "T"<<endl;
		  Ts+=1;
		  int T_to=species;
		  while (T_to==species) T_to=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(N); 	 
		  vector <long int> family;
		  long int father=population_of_genes[species][gene];
		  family.push_back(father);
		  long int daugther=next_gene;
		  next_gene++;
		  family.push_back(daugther);
		  long int son=next_gene;
		  next_gene++;
		  family.push_back(son);
		  vector < vector <long int> > fam_vec;
		  fam_vec.push_back(family);
		  gene_families.push_back(fam_vec);
		  gene_event++;

		  gene_event_times[father]=t;
		  gene_event_types[father]="T";

		  population_of_genes[species][gene]=daugther;
		  population_of_genes[T_to].push_back(son);
		  gene_count++;
		  //cout << "T."<<endl;

		}
	      else
		//L
		{
		  //cout << "L"<<endl;
		  Ls+=1;
		  int dj=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(population_of_genes[species].size());
		  population_of_genes[species].erase(population_of_genes[species].begin()+ dj );//XX is this what I wanted?
		  gene_count--;

		  //cout << "L."<<endl;

		}
	    }
	  else
	    //O
	    {
	      //cout << "O"<<endl;
	      Os+=1;
	      int species=RandomTools::giveIntRandomNumberBetweenZeroAndEntry(N);	    
	      population_of_genes[species].push_back(next_gene);
	      next_gene++;
	      gene_count++;

	      //cout << "O"<<endl;

	    }
	  
	}
      //species event happened
      else
	//S
	{
	  Ss+=1;
	  t=t_species_event;
	  //cout << "S"<<endl;

	  vector <long int> genes_in_daugther_species;
	  vector <long int> genes_in_son_species;
	  vector < vector <long int> > fam_vec;

	  for (vector <long int>::iterator gene=population_of_genes[births[species_event]].begin();gene!=population_of_genes[births[species_event]].end();gene++)
	    {
	      vector <long int> family;
	      long int father=(*gene);
	      family.push_back(father);
	      long int daugther=next_gene;
	      next_gene++;
	      family.push_back(daugther);
	      long int son=next_gene;
	      next_gene++;
	      family.push_back(son);
	      fam_vec.push_back(family);

	      genes_in_daugther_species.push_back(daugther);
	      genes_in_son_species.push_back(son);
	      gene_count++;
	      gene_event_types[father]="S";
	      gene_event_times[father]=t;
			  
	    }

	  gene_count-=population_of_genes[deaths[species_event]].size();

	  population_of_genes[births[species_event]].clear();
	  population_of_genes[births[species_event]]=genes_in_daugther_species;
	  population_of_genes[deaths[species_event]]=genes_in_son_species;

	  //cout << "S."<<endl;
	  

	  //XX oh my, what a waste of time ..
	  gene_families.push_back(fam_vec);
	  gene_event++;


      	  species_event++;
	  ++show_progress;
	}

    }

  cout << "#gene simulation ends "<< Ds << " Ds; "<< Ts << " Ts; " << Ls << " Ls; " << Os <<" Os; "<< Ss << " Ss."<<endl;


  //traceback gene stories
  map<long int,int> sampled_gene_counts;
  map<long int,string> gene_strings;
  map<long int,scalar_type> gene_age;
  int j=0;
  for (int i=0;i<n;i++) 
    {
      int sampled_i=sampled_population_indicies[i];
      long int extant_species = sampled_population[ i ];
      if (extant_species==population[ sampled_i ])
	for (vector<long int>::iterator git=population_of_genes[ sampled_i ].begin();git!=population_of_genes[ sampled_i ].end();git++)
	  {
	    long int extant_gene=(*git);
	    stringstream extant_gene_name;
	    extant_gene_name << i << "_" << j;// extant_species << "_" << extant_gene;
	    gene_strings[extant_gene]=extant_gene_name.str();
	    sampled_gene_counts[extant_gene]=1;  
	    gene_age[extant_gene]=0;
	    j++;
	  }
    }
  int gene_rank=0;

  cout << "#traceback begins.." <<endl;
  boost::progress_display show_trace_progress( gene_event);
  map <long int, string> suffix;
  for(vector <vector< vector<long int> > > ::reverse_iterator event_vec=gene_families.rbegin();event_vec!=gene_families.rend();event_vec++)
    {

      for(vector< vector<long int> >  ::iterator event=(*event_vec).begin();event!=(*event_vec).end();event++)
	{

	  long int father=(*event)[0];
	  long int daugther=(*event)[1];
	  long int son=(*event)[2];	      
	  scalar_type t_event=gene_event_times[father];

	  if (sampled_gene_counts[daugther]==1 and sampled_gene_counts[son]==1)
	    {
	      gene_rank++;
	      sampled_gene_counts[daugther]=0; 
	      sampled_gene_counts[son]=0; 
	      sampled_gene_counts[father]=1; 
	      stringstream sons_bl;
	      stringstream daugthers_bl;
	      stringstream gene_rank_bs;
	      gene_rank_bs << gene_event_types[father] << suffix[son] << suffix[daugther];
	      suffix[father]="";
	      sons_bl << t_event - gene_age[son];
	      daugthers_bl << t_event - gene_age[daugther];
	      gene_strings[father]="("+gene_strings[daugther]+":"+daugthers_bl.str()+","+gene_strings[son]+":"+sons_bl.str()+")"+gene_rank_bs.str();
	      gene_age[father]=t_event;
	      //lca=father;
	    }
	  else if (sampled_gene_counts[daugther]==1)
	    {
	      sampled_gene_counts[daugther]=0; 
	      sampled_gene_counts[father]=1; 
	      gene_strings[father]=gene_strings[daugther];
	      suffix[father]=suffix[daugther];
	      gene_age[father]=gene_age[daugther];

	    }
	  else if (sampled_gene_counts[son]==1)
	    {
	      sampled_gene_counts[son]=0; 
	      sampled_gene_counts[father]=1; 
	      gene_strings[father]=gene_strings[son];
	      suffix[father]=suffix[son];
	      if (gene_event_types[father]=="T") suffix[father]+="T";
	      gene_age[father]=gene_age[son];	  
	    }	
	}
      gene_event--;
      ++show_trace_progress;
    }
  
  int i=0;
  for (map<long int,int>::iterator git=sampled_gene_counts.begin();  git!=sampled_gene_counts.end(); git++)
    if ((*git).second==1 and gene_strings[(*git).first].find("(")!=string::npos )
      {
	stringstream fname;
	fname << "G_" << S_seed << "_" << G_seed << "_" << i << ".tree";
	ofstream gene_out( fname.str().c_str() );
	gene_out << gene_strings[(*git).first] << ";" << endl;
	i++;
      }

  
  return 1;
  

}
