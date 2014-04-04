#include "exODT.h"
#include <Bpp/Numeric/Random/RandomTools.h>


#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>


#include "exODT_sim.h"

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


int main(int args, char ** argv)
{

      BppApplication simulateSpAndGeneTrees(args, argv, "STRALE");
	simulateSpAndGeneTrees.startTimer();

  size_t N = ApplicationTools::getIntParameter("population.size",simulateSpAndGeneTrees.getParams(),1000);
  size_t n = ApplicationTools::getIntParameter("number.of.species",simulateSpAndGeneTrees.getParams(),10);
  size_t N_g = ApplicationTools::getIntParameter("number.of.genes",simulateSpAndGeneTrees.getParams(),100);

  double delta = ApplicationTools::getDoubleParameter("delta", simulateSpAndGeneTrees.getParams(), 0.01 );//0.20;
  double tau = ApplicationTools::getDoubleParameter("tau", simulateSpAndGeneTrees.getParams(), 0.01 );//0.31;
  double lambda = ApplicationTools::getDoubleParameter("lambda", simulateSpAndGeneTrees.getParams(), 0.1 );//0.39;

  	//If a seed is given
	int seed = ApplicationTools::getIntParameter("seed", simulateSpAndGeneTrees.getParams(), 0 );
	if (seed != 0 ) {
		RandomTools::setSeed ( seed ) ;
	}

  exODT_sim* simulation=new exODT_sim(N, seed);

//We then have to sample n species and get back a newick string of the represented phylogeny:

string Sstring=simulation->sample_species(n);

  stringstream fnameS;
  fnameS << "S" << ".tree";
  ofstream sp_out( fnameS.str().c_str() );
  sp_out << Sstring << ";" << endl;
  sp_out.close();

//We can after sampling ask for a vector of simulated gene tree newick strings

vector<string> Gstrings=simulation->simulate_gene_trees(2*N_g, delta, tau, lambda, 0, false, seed);

while (Gstrings.size() < N_g) {
  size_t todo = N_g - Gstrings.size() ;
  vector<string> Gstrings2=simulation->simulate_gene_trees(2*todo, delta, tau, lambda, 0, false, seed);
  VectorTools::append(Gstrings, Gstrings2);
}

  for (size_t i = 0; i<N_g ; ++i)
      {
	stringstream fname;
	fname << "G_" << i << ".tree";
	ofstream gene_out( fname.str().c_str() );
	gene_out << Gstrings[i] << ";" << endl;
	gene_out.close();
      }

  std::cout << "simulateSpAndGeneTrees's done. Bye." << std::endl;
  ApplicationTools::displayTime("Total execution time:");
  delete simulation;

//The simulate_gene_trees(..) function can be called as many times as needed. The represented species tree cannot be specified, but fixing the random seed will fix the species tree (both the represented and the unrepresented one). This is done like this:

//exODT_sim* simulation=new exODT_sim(N,S_seed);

//The parameters of the complete phylogeny simulation are the following:

/*exODT_sim* simulation=new exODT_sim(N,S_seed=-1,init_t=2);

    N: the number of species;
    S_seed: if =-1 it is taken from /dev/random.. ;
    init_t: the amount of time from the beginning of the simulation until the end in coalescent units.
*/
//The parameters of the gene tree simulation are:

/*vector<string> Gstrings=simulation->simulate_gene_trees(N_g,delta,tau,omega=0,only_root=false,G_seed=-1);

    N_g: number of gene families present at the start of the simulation;
    delta: D rate;
    tau: T rate;
    lambda: L rate;
    omega: uniform origination rate over the complete phylogeny;
    only_root: if true N_g gene families are present in the single species from which all n sampled species descend; if false N_g genes are present in each of the N species at the start of the simulation;
    G_Seed: if =-1 it is taken from /dev/random.. ;

The gene trees are time-like, i.e. are ultrametric and have branch lengths in coalescent units. The DTL events generating the family are not currently recorded. Doing this properly probably increases the complexity of the simulation. The current complexity is something like N^2 x (number of genes) in both time and space. I tested the a few things like extinction probability etc. that I could calculate analytically, so the simulation is probably correct.
*/

  
  return 1;
  

}
