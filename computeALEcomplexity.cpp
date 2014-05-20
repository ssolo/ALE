#include "exODT.h"
#include "ALE_util.h"
#include <boost/progress.hpp>

using namespace std;
using namespace bpp;

int main(int argc, char ** argv)
{
  //we need a species tree

  if (argc < 5) {
   std::cout << "\n\n\t\tUsage: computeALEcomplexity speciesTreeFile ALEfile delta tau lambda \n"<<std::endl;
       cout << "The program outputs statistics on an ale file as follows:\nale_file \t loglk \t NumberOfGenes\tNumberOfSpecies\tComputingTime1\tComputingTime2\tComputingTime3\tAverageComputingTime" << endl;
    std::cout.flush();
    exit(1);
  }
  
  string sname=argv[1];

  string Sstring;
  ifstream file_stream (sname.c_str());
  getline (file_stream,Sstring);

  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);

  exODT_model* model=new exODT_model();

  //XX
  //XX

  //exODT_sim* simulation=new exODT_sim(100,1010);
  
  scalar_type delta=atof(argv[3]);
  scalar_type tau=atof(argv[4]);
  scalar_type lambda=atof(argv[5]);

  //simulation->sample_species(10);

  
  //for (vector<string>::iterator it=simulation->gene_trees.begin();it!=simulation->gene_trees.end();it++)
  //cout << (*it) << endl;

  //cout << simulation->S_string << endl;

  model->set_model_parameter("min_D",1);
  model->set_model_parameter("grid_delta_t",0.01);
  model->set_model_parameter("DD",10);

  model->construct(Sstring);

  model->set_model_parameter("event_node",0);  
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("leaf_events",1);

  model->calculate_EGb();
  //Computing time:
  double ll;
  double time1, time2, time3;
  boost::timer * t = new boost::timer();

  ll = model->p(ale);
  time1 = t->elapsed();
  model->p(ale);
  time2 = t->elapsed();
  model->p(ale);
  time3 = t->elapsed();
  time3 = time3 - time2;
  time2 = time2 - time1;
  
  vector<string> leafNames = ale->getLeafNames();
  //Now, getting the number of species in the current ale_file
    std::map <std::string,scalar_type> scalar_parameter;
    std::map <std::string,std::string> string_parameter;
    string_parameter["gene_name_separators"]="_@";
    scalar_parameter["species_field"]=0;
    vector<string> speciesPresent;
  for (auto gene_name = leafNames.begin(); gene_name != leafNames.end(); ++gene_name) {
      vector <string> tokens;
      boost::split(tokens,*gene_name,boost::is_any_of(string_parameter["gene_name_separators"]),boost::token_compress_on);
      string species_name;
      if ((int)scalar_parameter["species_field"]==-1)
        species_name=tokens[tokens.size()-1];     
      else
        species_name=tokens[(int)scalar_parameter["species_field"]];      
      speciesPresent.push_back( species_name );
  }
  size_t numSpecies = (VectorTools::unique<string>(speciesPresent)).size();
      // gene<->species mapping

    cout << ale_file << "\t" << ll << "\t" << leafNames.size() << "\t" << numSpecies << "\t" << time1 << "\t" << time2 << "\t" << time3 << "\t" << (time1+time2+time3)/3 << endl;
    cout.flush();


  return 1;
  

}

