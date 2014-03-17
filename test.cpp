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

  model->set_model_parameter("min_D",20);
  model->set_model_parameter("grid_delta_t",0.005);
  model->set_model_parameter("DD",10);

  model->construct(Sstring);
  model->set_model_parameter("gene_name_separators", "-");

  model->set_model_parameter("event_node",0);  
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("leaf_events",1);

  model->calculate_EGb();
  //cout << model->p(ale) << endl;
  
  pair<string, scalar_type> res = model->p_MLRec(ale);    
  cout << res.first <<endl;
  cout << endl;
  cout << res.second <<endl;
  cout << endl;
  cout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
  cout <<"Total \t"<< model->MLRec_events["D"] << "\t" << model->MLRec_events["T"] << "\t" << model->MLRec_events["L"]<< "\t" << model->MLRec_events["S"] <<endl;      

  //for (int i=0;i<1;i++) cout << model->sample(false) << endl;
  return 1;
  

}

