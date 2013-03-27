#include "exODT.h"
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
  //we need an ale
  string ale_name=argv[2];
  approx_posterior * ale=load_ALE_from_file(ale_name);

  // initilaize the exODT model using some initial DTL rates
  exODT_model* model=new exODT_model();
  int D=atoi(argv[6]);
  model->set_model_parameter("D",D);
  int DD=atoi(argv[7]);
  model->set_model_parameter("DD",DD);
  model->construct(Sstring);
  model->set_model_parameter("event_node",1);
  scalar_type delta=atof(argv[3]);
  scalar_type tau=atof(argv[4]);
  scalar_type lambda=atof(argv[5]);
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau",tau);
  model->set_model_parameter("lambda",lambda);

  model->set_model_parameter("event_node",0);
  model->calculate_EGb();

  pair<string, scalar_type> res =  model->p_MLRec(ale); 
  cout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
  cout << res.first << endl;
  cout <<"Total \t"<< model->MLRec_events["D"] << "\t" << model->MLRec_events["T"] << "\t" << model->MLRec_events["L"]<< "\t" << model->MLRec_events["S"] <<endl;    

  return 1;
  //model->set_model_parameter("event_node",1);

  for (tau=delta;tau<1e4;tau*=sqrt(10))
    {
      model->set_model_parameter("delta",tau);
      // calculate single gene propagtaion and extinction functions
      model->calculate_EGb();
      cout << tau << " " << model->p(ale) << endl;
    }
  

}

