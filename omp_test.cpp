#include "exODT.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
int main(int argc, char ** argv)
{
  //we need a species tree
  string Sstring;
  ifstream file_stream ("example_data/cy36_green.tree");
  getline (file_stream,Sstring);
  //we need an ale
  string ale_name=argv[1];
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
  cout << model->p(ale) << endl;
  // the openMP part should be added in model.cpp in exODT_model::p( approx_posterior *)
  // I guess you have to add the openMP bit into model.cpp while making sure the above number does not change 
  // and also make sure that the tables produed during the calculations work ..
  RandomTools::setSeed(20110426);
  cout << model->sample() << endl;
  // this can be done with the stochastic backtrace, since we fixed the random seed, the above reconciled tree should alos not change..  
  cout << model->p(ale) << endl;
  

}

