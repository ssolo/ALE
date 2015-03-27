
#include "exODT.h"
#include "ALE_util.h"

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>

using namespace std;
using namespace bpp;

int main(int argc, char ** argv)
{


  //we need a dated species tree in newick format
  string Sstring;
  ifstream file_stream_S (argv[1]);
  getline (file_stream_S,Sstring);
  cout << "Read species tree from: " << argv[1] <<".."<<endl;
  //we need an .ale file containing observed conditional clade probabilities
  //cf. ALEobserve
  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);
  cout << "Read summary of tree sample for "<<ale->observations<<" trees from: " << ale_file <<".."<<endl;

  //we initialise a coarse grained reconciliation model for calculating the sum
  exODT_model* model=new exODT_model();

  if (argc>3)
    model->set_model_parameter("gene_name_separators", argv[3]);
  model->set_model_parameter("BOOT_STRAP_LABLES","yes");
  model->construct_undated(Sstring);


  //a set of inital rates
  scalar_type delta=0.01,tau=0.01,lambda=0.1;  
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->calculate_undatedEs();

  tree_type * T1=TreeTemplateTools::parenthesisToTree(ale->constructor_string,false);
  vector <Node *> nodes1=T1->getLeaves();
  map<string,int> names;
  for (vector <Node *>::iterator it=nodes1.begin();it!=nodes1.end();it++)
    {
      vector <string> tokens;
      string name=(*it)->getName();
      boost::split(tokens,name,boost::is_any_of("_"),boost::token_compress_on);
      names[tokens[0]]+=1;
    }
  

  boost::timer * t = new boost::timer();
  string outname=ale_file+".utimes"; 
  ofstream fout( outname.c_str() );
  model->pun(ale);
  fout << t->elapsed() << "\t";
  fout << ale->Dip_counts.size() << "\t";
  fout <<  T1->getNumberOfLeaves() << "\t";
  fout <<  names.size() << "\t";      
  fout << ale_file << endl;
  
  return 0;
}

