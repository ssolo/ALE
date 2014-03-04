#include "exODT.h"
#include "exODT_sim.h"

#include "ALE_util.h"
#include <omp.h>

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  //we need a species tree

  string sname=argv[1];

  string Sstring;
  ifstream file_stream (sname.c_str());
  getline (file_stream,Sstring);

  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);

  exODT_model* model=new exODT_model();
  exODT_sim* simulation=new exODT_sim(100,1010);
  
  simulation->sample_species(10);

  simulation->simulate_gene_trees(10,0.05,0.05,0.05);
  
  for (vector<string>::iterator it=simulation->gene_trees.begin();it!=simulation->gene_trees.end();it++)
    cout << (*it) << endl;

  cout << simulation->S_string << endl;

  model->set_model_parameter("min_D",1);
  model->set_model_parameter("grid_delta_t",0.5);
  model->set_model_parameter("DD",30);

  model->construct(Sstring);

  model->set_model_parameter("event_node",0);  
  model->set_model_parameter("delta", atof(argv[3]));
  model->set_model_parameter("tau", atof(argv[4]));
  model->set_model_parameter("lambda", atof(argv[5]));

  model->calculate_EGb();
  cout << model->p(ale) << endl;

  return 1;
  

}

