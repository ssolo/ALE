#include "exODT.h"
#include "ALE_util.h"
#include <omp.h>

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  cout << "changing OMP_NUM_THREADS should leave the results unchanged!" << endl;
  //we need a species tree
  string Sstring;
  ifstream file_stream ("example_data/cy36_green.tree");
  getline (file_stream,Sstring);
  //we need an ale
  string ale_name="example_data/sc_cy36HBG285662.ale";
  if (argc>1) ale_name=argv[1];
  approx_posterior * ale=load_ALE_from_file(ale_name);

  // initilaize the exODT model using some initial DTL rates
  exODT_model* model=new exODT_model();
  model->set_model_parameter("D",3);
  model->set_model_parameter("DD",10);
  model->construct(Sstring);
  model->set_model_parameter("event_node",0);
  model->set_model_parameter("delta",0.1);
  model->set_model_parameter("tau",0.1);
  model->set_model_parameter("lambda",0.2);
  // calculate single gene propagtaion and extinction functions
  model->calculate_EGb();
  // calculate joint ALE*exODT likelihood summed over all reconcilation and all tree toplopgies

  // the openMP part should be added in model.cpp in exODT_model::p( approx_posterior *)
  // I guess you have to add the openMP bit into model.cpp while making sure the above number does not change 
  // and also make sure that the tables produed during the calculations work ..
  // this can be done with the stochastic backtrace, since we fixed the random seed, the above reconciled tree should alos not change..  
  // calling p(ale) after sample when currrently unactivated pragma with label //p3 is active crashes 

  scalar_type t_0,t_1;
  omp_set_num_threads(1);

  cout << endl << "trying OMP_NUM_THREADS=1" << endl; 
  t_0=omp_get_wtime();
  cout << model->p(ale) << endl;
  t_1=omp_get_wtime();
  cout << t_1-t_0 <<"s"<< endl << endl;

  RandomTools::setSeed(20110426);
  cout << model->sample() << endl;

  omp_set_num_threads(2);
  cout << endl << "trying OMP_NUM_THREADS=2" << endl; 
  t_0=omp_get_wtime();
  cout << model->p(ale) << endl;
  t_1=omp_get_wtime();
  cout << t_1-t_0 <<"s"<< endl << endl;

  RandomTools::setSeed(20110426);
  cout << model->sample() << endl;

  omp_set_num_threads(4);
  cout << endl << "trying OMP_NUM_THREADS=4" << endl; 
  t_0=omp_get_wtime();
  cout << model->p(ale) << endl;
  t_1=omp_get_wtime();
  cout << t_1-t_0 <<"s"<< endl << endl;

  RandomTools::setSeed(20110426);
  cout << model->sample() << endl;

  omp_set_num_threads(8);
  cout << endl << "trying OMP_NUM_THREADS=8" << endl; 
  t_0=omp_get_wtime();
  cout << model->p(ale) << endl;
  t_1=omp_get_wtime();
  cout << t_1-t_0 <<"s"<< endl << endl;

  RandomTools::setSeed(20110426);
  cout << model->sample() << endl;

  omp_set_num_threads(12);
  cout << endl << "trying OMP_NUM_THREADS=12" << endl;
  t_0=omp_get_wtime();
  cout << model->p(ale) << endl;
  t_1=omp_get_wtime();
  cout << t_1-t_0 <<"s"<< endl << endl;

  RandomTools::setSeed(20110426);
  cout << model->sample() << endl;
  

}

