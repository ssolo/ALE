
#include "ALE_util.h"
#include "exODT.h"
#include <boost/progress.hpp>

#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>

using namespace std;
using namespace bpp;

class p_fun : public virtual Function, public AbstractParametrizable {
private:
  double fval_;
  bool delta_fixed;
  bool tau_fixed;
  bool lambda_fixed;
  exODT_model *model_pointer;
  approx_posterior *ale_pointer;

public:
  p_fun(exODT_model *model, approx_posterior *ale, double delta_start = 0.01,
        double tau_start = 0.01,
        double lambda_start = 0.1 //,double sigma_hat_start=1.
        ,
        bool delta_fixed_in = false, bool tau_fixed_in = false,
        bool lambda_fixed_in = false)
      : AbstractParametrizable(""), fval_(0), model_pointer(model),
        ale_pointer(ale) {
    // We declare parameters here:
    //   IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
    IntervalConstraint *constraint =
        new IntervalConstraint(1e-5, 1, true, true);
    delta_fixed = delta_fixed_in;
    tau_fixed = tau_fixed_in;
    lambda_fixed = lambda_fixed_in;

    if (not delta_fixed) {
      addParameter_(new Parameter("delta", delta_start, constraint));
      cout << "#optimizing delta rate" << endl;
    }
    if (not tau_fixed) {
      addParameter_(new Parameter("tau", tau_start, constraint));
      cout << "#optimizing tau rate" << endl;
    }
    if (not lambda_fixed) {
      addParameter_(new Parameter("lambda", lambda_start, constraint));
      cout << "#optimizing lambda rate" << endl;
    }
  }

  p_fun *clone() const { return new p_fun(*this); }

public:
  void setParameters(const ParameterList &pl) throw(ParameterNotFoundException,
                                                    ConstraintException,
                                                    Exception) {
    matchParametersValues(pl);
  }
  double getValue() const throw(Exception) { return fval_; }
  void fireParameterChanged(const ParameterList &pl) {
    double delta, tau, lambda;
    if (not delta_fixed) {
      delta = getParameterValue("delta");
      model_pointer->set_model_parameter("delta", delta);
    }
    if (not tau_fixed) {
      tau = getParameterValue("tau");
      model_pointer->set_model_parameter("tau", tau);
    }
    if (not lambda_fixed) {
      lambda = getParameterValue("lambda");
      model_pointer->set_model_parameter("lambda", lambda);
    }

    model_pointer->calculate_EGb();
    double y = -log(model_pointer->p(ale_pointer));
    cout << "delta=" << delta << "\t tau=" << tau << "\t lambda=" << lambda
         << "\t ll=" << -y << endl;
    fval_ = y;
  }
};

int main(int argc, char **argv) {
  cout << "ALEml using ALE v" << ALE_VERSION << endl;

  if (argc < 3) {
    cout << "\nUsage:\n ./ALEml species_tree.newick gene_tree_sample.ale "
            "sample=number_of_samples separators=gene_name_separator "
            "O_R=OriginationAtRoot delta=DuplicationRate tau=TransferRate "
            "lambda=LossRate beta=weight_of_sequence_evidence "
            "fraction_missing=file_with_fraction_of_missing_genes_per_species "
            "output_species_tree=n S_branch_lengths:root_length "
            "rate_mutiplier:rate_name:branch_id:value"
         << endl;
    cout << "\n1st example: we fix the DTL values and do not perform any "
            "optimization \n ./ALEml species_tree.newick gene_tree_sample.ale "
            "sample=100 separators=_ delta=0.05 tau=0.1 lambda=0.2 "
         << endl;
    cout << "\n2nd example: we fix the T value to 0 to get a DL-only model and "
            "optimize the DL parameters \n ./ALEml species_tree.newick "
            "gene_tree_sample.ale sample=100 separators=_ tau=0\n"
         << endl;

    return 0;
  }

  // we need a dared species tree in newick format
  string Sstring;
  string S_treefile = argv[1];
  if (!fexists(argv[1])) {
    cout << "Error, file " << argv[1] << " does not seem accessible." << endl;
    exit(1);
  }
  ifstream file_stream_S(argv[1]);
  getline(file_stream_S, Sstring);
  cout << "Read species tree from: " << argv[1] << ".." << endl;
  // we need an .ale file containing observed conditional clade probabilities
  // cf. ALEobserve
  string ale_file = argv[2];
  approx_posterior *ale;
  ale = load_ALE_from_file(ale_file);
  cout << "Read summary of tree sample for " << ale->observations
       << " trees from: " << ale_file << ".." << endl;
  // Getting the radical for output files:
  vector<string> tokens;
  boost::split(tokens, ale_file, boost::is_any_of("/"),
               boost::token_compress_on);
  ale_file = tokens[tokens.size() - 1];
  boost::split(tokens, S_treefile, boost::is_any_of("/"),
               boost::token_compress_on);
  ale_file = tokens[tokens.size() - 1] + "_" + ale_file;

  // we initialise a coarse grained reconciliation model for calculating the sum
  exODT_model *model = new exODT_model();

  scalar_type samples = 100;
  scalar_type O_R = 1, beta = 1;
  bool delta_fixed = false;
  bool tau_fixed = false;
  bool lambda_fixed = false;
  scalar_type delta = 1e-2, tau = 1e-2, lambda = 1e-1;
  string fractionMissingFile = "";
  string outputSpeciesTree = "";
  // ############################################################
  for (int i = 3; i < argc; i++) {
    string next_field = argv[i];

    vector<string> tokens;
    boost::split(tokens, next_field, boost::is_any_of("=:"),
                 boost::token_compress_on);
    if (tokens[0] == "sample")
      samples = atoi(tokens[1].c_str());
    else if (tokens[0] == "separators")
      model->set_model_parameter("gene_name_separators", tokens[1]);
    else if (tokens[0] == "delta") {
      delta = atof(tokens[1].c_str());
      delta_fixed = true;
      cout << "# delta fixed to " << delta << endl;
    } else if (tokens[0] == "tau") {
      tau = atof(tokens[1].c_str());
      tau_fixed = true;
      cout << "# tau fixed to " << tau << endl;
    } else if (tokens[0] == "lambda") {
      lambda = atof(tokens[1].c_str());
      lambda_fixed = true;
      cout << "# lambda fixed to " << lambda << endl;
    } else if (tokens[0] == "O_R") {
      O_R = atof(tokens[1].c_str());
      cout << "# NOT YET IMPLEMENTED O_R set to " << O_R << endl;
    } else if (tokens[0] == "beta") {
      beta = atof(tokens[1].c_str());
      cout << "# beta set to " << beta << endl;
    } else if (tokens[0] == "fraction_missing") {
      fractionMissingFile = tokens[1];
      cout << "# NOT YET IMPLEMENTED File containing fractions of missing "
              "genes set to "
           << fractionMissingFile << endl;
    } else if (tokens[0] == "output_species_tree") {
      std::string valArg = boost::algorithm::to_lower_copy(tokens[1]);
      if (valArg == "y" || valArg == "ye" || valArg == "yes") {
        outputSpeciesTree = ale_file + ".spTree";
        cout << "# outputting the annotated species tree to "
             << outputSpeciesTree << endl;
      } else {
        cout << "# NOT outputting the annotated species tree to "
             << outputSpeciesTree << endl;
      }
    }
  }
  // ############################################################

  int D = 2;
  model->set_model_parameter("BOOTSTRAP_LABELS", "yes");

  model->set_model_parameter("min_D", D);
  model->set_model_parameter("grid_delta_t", 0.1);
  model->construct(Sstring);

  model->set_model_parameter("event_node", 0);
  model->set_model_parameter("leaf_events", 1);
  model->set_model_parameter("N", 1);

  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);
  model->set_model_parameter("sigma_hat", 1);

  // calculate_EGb() must always be called after changing rates to calculate E-s
  // and G-s cf. http://arxiv.org/abs/1211.4606
  model->calculate_EGb();

  cout << "Reconciliation model initialised, starting DTL rate optimisation"
       << ".." << endl;

  // we use the Nelder–Mead method implemented in Bio++
  Function *f = new p_fun(model, ale, delta, tau, lambda, delta_fixed,
                          tau_fixed, lambda_fixed);
  Optimizer *optimizer = new DownhillSimplexMethod(f);

  optimizer->setProfiler(0);
  optimizer->setMessageHandler(0);
  optimizer->setVerbose(2);

  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  optimizer->init(f->getParameters()); // Here we optimize all parameters, and
                                       // start with the default values.
  scalar_type mlll;
  if (not(delta_fixed and tau_fixed and lambda_fixed)) {
    cout << "#optimizing rates" << endl;
    optimizer->optimize();
    if (not delta_fixed)
      delta = optimizer->getParameterValue("delta");
    if (not tau_fixed)
      tau = optimizer->getParameterValue("tau");
    if (not lambda_fixed)
      lambda = optimizer->getParameterValue("lambda");
    mlll = -optimizer->getFunctionValue();
  } else {
    mlll = log(model->p(ale));
  }
  cout << endl
       << "ML rates: "
       << " delta=" << delta << "; tau=" << tau
       << "; lambda=" << lambda //<<"; sigma="<<sigma_hat
       << "." << endl;

  cout << "LL=" << mlll << endl;

  cout << "Sampling reconciled gene trees.." << endl;
  vector<string> sample_strings;
  vector<Tree *> sample_trees;
  boost::progress_display pd(samples);

  for (int i = 0; i < samples; i++) {
    ++pd;
    string sample_tree = model->sample(false);
    sample_strings.push_back(sample_tree);

    if (ale->last_leafset_id > 3) {

      tree_type *G = TreeTemplateTools::parenthesisToTree(sample_tree, false);
      vector<Node *> leaves = G->getLeaves();
      for (vector<Node *>::iterator it = leaves.begin(); it != leaves.end();
           it++) {
        string name = (*it)->getName();
        vector<string> tokens;
        boost::split(tokens, name, boost::is_any_of(".@"),
                     boost::token_compress_on);
        (*it)->setName(tokens[0]);
        tokens.clear();
      }
      leaves.clear();
      sample_trees.push_back(G);
    }
  }
  /*cout << "Calculating ML reconciled gene tree.."<<endl;
  pair<string, scalar_type> res = model->p_MLRec(ale);
  //and output it..
  */
  string outname = ale_file + ".ml_rec";
  ofstream fout(outname.c_str());
  fout << "#ALEml using ALE v" << ALE_VERSION
       << " by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;" << endl
       << endl;
  fout << "S:\t" << model->string_parameter["S_with_ranks"] << endl;
  fout << endl;
  fout << "Input ale from:\t" << ale_file << endl;
  fout << ">logl: " << mlll << endl;
  fout << "rate of\t Duplications\tTransfers\tLosses" << endl;
  fout << "ML \t" << delta << "\t" << tau << "\t"
       << lambda //<< "'t" << sigma_hat
       << endl;
  fout << endl;
  fout << samples << " reconciled G-s:\n" << endl;
  for (int i = 0; i < samples; i++) {
    fout << sample_strings[i] << endl;
  }

  // fout << "reconciled G:\t"<< res.first <<endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" << endl;
  fout << "Total \t" << model->MLRec_events["D"] / samples << "\t"
       << model->MLRec_events["T"] / samples << "\t"
       << model->MLRec_events["L"] / samples << "\t"
       << model->MLRec_events["S"] / samples << endl;
  fout << endl;
  fout << "# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies"
       << endl;
  fout << model->counts_string(samples);
  fout.close();

  // Outputting the species tree to its own file:
  if (outputSpeciesTree != "") {
    ofstream fout(outputSpeciesTree.c_str());
    fout << model->string_parameter["S_with_ranks"] << endl;
    fout.close();
  }

  cout << "Results in: " << outname << endl;
  if (ale->last_leafset_id > 3) {
    cout << "Calculating MRP consensus tree." << endl;
    Tree *con_tree = TreeTools::MRP(sample_trees);

    string con_name = ale_file + ".cons_tree";

    ofstream con_out(con_name.c_str());
    con_out << "#ALEsample using ALE v" << ALE_VERSION
            << " by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;" << endl;
    TreeTools::computeBootstrapValues(*con_tree, sample_trees);
    string con_tree_sup = TreeTemplateTools::treeToParenthesis(*con_tree);
    con_out << con_tree_sup << endl;
    cout << endl << "Consensus tree in " << con_name << endl;
  }

  string t_name = ale_file + ".Ts";
  ofstream tout(t_name.c_str());
  tout << "#Transfer & duplications tokens: \n";
  tout << "#D|rank|named_branch|g_id\n;";
  tout << "#T(|rank|t|named_branch|g_id)>(|rank|t|named_branch|gp_id)>.. or \n";
  tout << "#>(|rank|t|named_branch|-1), where g_id=-1 is the root of G.\n";

  for (std::vector<std::string>::iterator it = model->Ttokens.begin();
       it != model->Ttokens.end(); it++)
    tout << (*it) << endl;
  tout.close();

  cout << "Transfers in: " << t_name << endl;
  return 0;
}
