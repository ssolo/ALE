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


  cout << "Read species tree from: " << argv[1] <<".."<<endl;
  //we need an .ale file containing observed conditional clade probabilities
  //cf. ALEobserve
  cout << "Read summary of tree sample for "<<ale->observations<<" trees from: " << ale_name <<".. with :" << ale->count_trees() << " possible amalgamations .." << endl << endl ;

  cout << "and the most liely tree is.."<< endl;

  // initilaize the exODT model using some initial DTL rates
  exODT_model* model=new exODT_model();
  int D=1; // this is the simplest parsimony like setting of one event node per slice.. 
  model->set_model_parameter("D",D);

  model->construct(Sstring);

  scalar_type delta=0.01;
  scalar_type tau=0.01;
  scalar_type lambda=0.01;
  model->set_model_parameter("delta",delta);
  model->set_model_parameter("tau",tau);
  model->set_model_parameter("lambda",lambda);

  // likelihood only calculation, the function E and G in the paper.. 
  model->calculate_EGb();
  
  //find ML reconciled gene tree:
  pair<string, scalar_type> res =  model->p_MLRec(ale); 
  //output tree
  cout << res.first << endl;
  //output number of events
  cout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl; 
  cout <<"Total \t"<< model->MLRec_events["D"] << "\t" << model->MLRec_events["T"] << "\t" << model->MLRec_events["L"]<< "\t" << model->MLRec_events["S"] <<endl;    

  //ususally one needs to run instead of model->p_MLRec(ale) the below: 

  //model->calculate_EGb();
  //cout << model->p(ale) << endl;

  //so maybe makes more sense to compare to that .. 

  return 1;
  

}

