//all code by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
#include "ALE.h"
#include <Bpp/Numeric/Random/RandomTools.h>

class exODT_sim
{
 public:
  long int S_seed,G_seed;
  int N,n;

  scalar_type init_t,sigma,lca_age;

  std::vector<std::string> gene_trees;

  std::string S_string,R_string;

  exODT_sim( int N_in, long int S_seed_in=-1, scalar_type init_t_in=2, scalar_type sigma=-1);
    
  std::string sample_species(int n_in);

  std::vector<std::string> simulate_gene_trees(int G_n,scalar_type delta,scalar_type tau,scalar_type lambda,scalar_type omega=0,bool only_root=true, long int G_seed_in=-1,bool event_string=false);
    
  ~exODT_sim()
    {
      population.clear();

      for (std::vector<std::vector<long int> >::iterator it=families.begin();it!=families.end();it++)
	(*it).clear();
      families.clear();

      event_times.clear();

      births.clear();
      
      deaths.clear();

      sampled_population.clear();

      sampled_population_indicies.clear();
      
    }


 private:
  std::vector <long int> population;//del-loc
  long int next_index;
  long long species_event;
  long long number_of_species_events;

  std::vector<std::vector<long int> > families;//del-loc
  std::map <long long,scalar_type > event_times;//del-loc

  std::vector < int > births;//del-loc
  std::vector < int > deaths;//del-loc

  std::vector<long int> sampled_population;//del-loc
  std::vector <int> sampled_population_indicies;//del-loc  

  long int lca;
  int fca;


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
  }
  
};
