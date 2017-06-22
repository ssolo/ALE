#include "ALE.h"
#include "ALE_util.h"
#include "exODT.h"

using namespace std;
using namespace bpp;


string readTreeFromFile(string fname) {
  if (!fexists(fname)) {
    cout << "Error, file "<<fname << " does not seem accessible." << endl;
    exit(1);
  }
  ifstream file_stream (fname.c_str());
  string tree_str="";
  if (file_stream.is_open())  //  ########## read trees ############
  {
    while (! file_stream.eof())
    {
      string line;
      getline (file_stream,line);
      if (line.find("(")!=line.npos )
      {
        tree_str = line;
      }
    }
  }
  return tree_str;
}

int main(int argc, char ** argv)
{
  cout << "ALEestimate using ALE v"<< ALE_VERSION <<endl;

  if (argc<2)
  {
    cout << "usage:\n ./ALEestimate species_tree_file gene_tree_file separators=gene_name_separator O_R=OriginationAtRoot delta=DuplicationRate tau=TransferRate lambda=LossRate beta=weight_of_sequence_evidence outputFiles=n" << endl;
    return 1;
  }

  // Species tree business
  string species_tree_file=argv[1];
  boost::trim(species_tree_file);
  //Reading the species tree from within the file
  string species_tree_str = readTreeFromFile( species_tree_file );
  cout << "\n\tRead species tree from: " << argv[1] <<endl;

  // Gene tree business
  string gene_tree_file=argv[2];
  boost::trim(gene_tree_file);
  string head=gene_tree_file;
  string ale_name=head+".ale";
  //Reading the gene tree from within the file
  string gene_tree_str = readTreeFromFile( gene_tree_file );
  approx_posterior* ale = new approx_posterior(gene_tree_str);
  vector<string> gene_tree_strs ; // Silly: we need to produce a vector with a single element...
  gene_tree_strs.push_back(gene_tree_str);
  ale->observation(gene_tree_strs);
  cout << "\n\tObserved "<< ale->observations << " gene tree(s) from: " <<  argv[2] <<endl ;

  //we initialise a coarse grained reconciliation model for calculating the sum
  exODT_model* model=new exODT_model();

  // Getting the other options
  scalar_type samples=100;
  scalar_type O_R=1,beta=1;
  scalar_type delta=0.01,tau=0.01,lambda=0.1;
  string fractionMissingFile = "";
  bool outputyn = false;
  for (int i=3;i<argc;i++)
  {
    string next_field=argv[i];
    vector <string> tokens;
    boost::split(tokens,next_field,boost::is_any_of("="),boost::token_compress_on);
    if (tokens[0]=="sample")
    samples=atoi(tokens[1].c_str());
    else if (tokens[0]=="separators")
    model->set_model_parameter("gene_name_separators", tokens[1]);
    else if (tokens[0]=="delta")
    {
      delta=atof(tokens[1].c_str());
      cout << "\n\tDelta fixed to " << delta << endl;
    }
    else if (tokens[0]=="tau")
    {
      tau=atof(tokens[1].c_str());
      cout << "\n\tTau fixed to " << tau << endl;
    }
    else if (tokens[0]=="lambda")
    {
      lambda=atof(tokens[1].c_str());
      cout << "Lambda fixed to " << lambda << endl;
    }
    else if (tokens[0]=="O_R")
    {
      O_R=atof(tokens[1].c_str());
      cout << "\n\tO_R set to " << O_R << endl;
    }
    else if (tokens[0]=="beta")
    {
      beta=atof(tokens[1].c_str());
      cout << "\n\tBeta set to " << beta << endl;
    }
    else if (tokens[0]=="fraction_missing")
    {
      fractionMissingFile=tokens[1];
      cout << "\n\tFile containing fractions of missing genes set to " << fractionMissingFile << endl;
    }
    else if (tokens[0]=="outputFiles")
    {
      if (tokens[1] == "y" || tokens[1] == "yes" || tokens[1] == "Y" || tokens[1] == "YES") {
        outputyn = true;
      }
    }
  }

  // Constructing the ALE_undated object and computing the logLk.

  model->set_model_parameter("BOOTSTRAP_LABELS","yes");
  model->construct_undated(species_tree_str, fractionMissingFile);

  model->set_model_parameter("seq_beta", beta);
  model->set_model_parameter("O_R", O_R);
  //a set of inital rates
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);

  //calculate_EGb() must always be called after changing rates to calculate E-s and G-s
  //cf. http://arxiv.org/abs/1211.4606
  model->calculate_undatedEs();
  double loglk = log(model->pun(ale, true));
  cout << "\n\tReconciliation model likelihood computed, logLk: " <<loglk<<endl;

  // Output
  if (outputyn) {
    cout << "\n\tSampling reconciled gene trees.."<<endl;
    vector <string> sample_strings;
    vector <Tree*> sample_trees;
    boost::progress_display pd( samples );

    for (int i=0;i<samples;i++)
    {
      ++pd;
      string sample_tree=model->sample_undated();
      sample_strings.push_back(sample_tree);
      if (ale->last_leafset_id>3)
      {

        tree_type * G=TreeTemplateTools::parenthesisToTree(sample_tree,false);

        vector<Node*> leaves = G->getLeaves();
        for (vector<Node*>::iterator it=leaves.begin();it!=leaves.end();it++ )
        {
          string name=(*it)->getName();
          vector<string> tokens;
          boost::split(tokens,name,boost::is_any_of(".@"),boost::token_compress_on);
          (*it)->setName(tokens[0]);
          tokens.clear();
        }
        leaves.clear();
        sample_trees.push_back(G);
      }
    }
    vector<string> tokens;
    boost::split(tokens,gene_tree_file,boost::is_any_of("/"),boost::token_compress_on);
    ale_name=tokens[tokens.size()-1];
    string outname=ale_name+".uml_rec";
    ofstream fout( outname.c_str() );
    fout <<  "#ALEevaluate using ALE v"<< ALE_VERSION <<"; CC BY-SA 3.0;"<<endl<<endl;
    fout << "S:\t"<<model->string_parameter["S_with_ranks"] <<endl;
    fout << endl;
    fout << "Gene tree from:\t"<<gene_tree_file<<endl;
    fout << ">logl: " << loglk << endl;
    fout << "rate of\t Duplications\tTransfers\tLosses" <<endl;
    fout << "\t"<< delta << "\t" << tau << "\t" << lambda //<< "'t" << sigma_hat
    << endl;
    fout << endl;
    fout << samples << " reconciled G-s:\n"<<endl;
    for (int i=0;i<samples;i++)
    {
      fout<<sample_strings[i]<<endl;
    }

    fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl;
    fout <<"Total \t"<< model->MLRec_events["D"]/samples << "\t" << model->MLRec_events["T"]/samples << "\t" << model->MLRec_events["L"]/samples<< "\t" << model->MLRec_events["S"]/samples <<endl;
    fout << endl;
    fout << "# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies" <<endl;
    fout << model->counts_string_undated(samples);
    fout.close();

    cout << "Results in: " << outname << endl;
    if (ale->last_leafset_id>3)
    {
      cout << "Calculating consensus tree."<<endl;
      Tree* con_tree= TreeTools::thresholdConsensus(sample_trees,0.5);

      string con_name=ale_name+".ucons_tree";

      ofstream con_out( con_name.c_str() );
      con_out <<  "#ALEsample using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl;
      TreeTools::computeBootstrapValues(*con_tree,sample_trees);
      string con_tree_sup=TreeTemplateTools::treeToParenthesis(*con_tree);
      con_out << con_tree_sup << endl;
      con_out.close();
      cout << endl<< "Consensus tree in " << con_name<< endl;
    }

    string t_name=ale_name+".uTs";
    ofstream tout( t_name.c_str() );
    tout <<"#from\tto\tfreq.\n";

    for (int e=0;e<model->last_branch;e++)
    for (int f=0;f<model->last_branch;f++)
    if  (model->T_to_from[e][f]>0)
    {
      if (e<model->last_leaf)
      tout << "\t" << model->node_name[model->id_nodes[e]];
      else
      tout << "\t" << e;
      if (f<model->last_leaf)
      tout << "\t" << model->node_name[model->id_nodes[f]];
      else
      tout << "\t" << f;
      tout << "\t" << model->T_to_from[e][f]/samples <<  endl;
    }
    tout.close();
    cout << "Transfers in: " << t_name << endl;
  }

  return 0;
}
