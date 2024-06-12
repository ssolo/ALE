#include "ALE_util.h"
#include "exODT.h"

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
  bool DT_fixed;
  bool MLOR;

  bool no_T = false;
  exODT_model *model_pointer;
  approx_posterior *ale_pointer;

public:
  p_fun(exODT_model *model, approx_posterior *ale, vector<int> ml_branch_ids,
        vector<string> ml_ratetype_names, double delta_start = 0.1,
        double tau_start = 0.1,
        double lambda_start = 0.5 //,double sigma_hat_start=1.
        ,
        bool delta_fixed_in = false, bool tau_fixed_in = false,
        bool lambda_fixed_in = false, bool DT_fixed_in = false,
        bool MLOR_in = false)
      : AbstractParametrizable(""), fval_(0), model_pointer(model),
        ale_pointer(ale) {
    // We declare parameters here:
    //    IncludingInterval* constraint = new IncludingInterval(1e-6, 10-1e-6);
    IntervalConstraint *constraint;
    if (no_T)
      constraint = new IntervalConstraint(1e-6, 100 - 1e-7, true, true);
    else
      constraint = new IntervalConstraint(1e-10, 100 - 1e-7, true, true);

    IntervalConstraint *rate_multiplier_constraint =
        new IntervalConstraint(1e-7, 10000 - 1e-7, true, true);
    delta_fixed = delta_fixed_in;
    tau_fixed = tau_fixed_in;
    lambda_fixed = lambda_fixed_in;
    DT_fixed = DT_fixed_in;
    MLOR = MLOR_in;
    if (tau_start < 1e-10)
      no_T = true;
    if (not delta_fixed and not DT_fixed) {
      addParameter_(new Parameter("delta", delta_start, constraint));
      cout << "#optimizing delta rate" << endl;
    }
    if (not tau_fixed and not DT_fixed) {
      addParameter_(new Parameter("tau", tau_start, constraint));
      cout << "#optimizing tau rate" << endl;
    }
    if (not lambda_fixed) {
      addParameter_(new Parameter("lambda", lambda_start, constraint));
      cout << "#optimizing lambda rate" << endl;
    }
    if (DT_fixed) {
      addParameter_(new Parameter("tau", tau_start, constraint));
      cout << "#optimizing delta and tau rates with fixed D/T ratio" << endl;
    }
    if (MLOR) {
      IntervalConstraint *OR_constraint =
          new IntervalConstraint(1e-10, 1000 - 1e-7, true, true);
      addParameter_(new Parameter("O_R", 1., OR_constraint));
      cout << "#optimizing O_R" << endl;
    }

    // vector<int> ml_branch_ids;
    // vector<string> ml_ratetype_names;

    for (int i = 0; i < ml_branch_ids.size(); i++) {
      int e = ml_branch_ids[i];
      public_ml_branch_ids.push_back(e);
      string name = ml_ratetype_names[i];
      public_ml_ratetype_names.push_back(name);
      stringstream branch;
      branch << e;
      addParameter_(new Parameter("rm_" + name + "_" + branch.str(), 1,
                                  rate_multiplier_constraint));
      cout << "#optimizing for branch " << e << " ratemultiplier " << name
           << endl;
    }
  }

  p_fun *clone() const { return new p_fun(*this); }

public:
  vector<int> public_ml_branch_ids;
  vector<string> public_ml_ratetype_names;

  void setParameters(const ParameterList &pl) throw(ParameterNotFoundException,
                                                    ConstraintException,
                                                    Exception) {
    matchParametersValues(pl);
  }
  double getValue() const throw(Exception) { return fval_; }
  void fireParameterChanged(const ParameterList &pl) {
    if (not delta_fixed and not DT_fixed) {
      double delta = getParameterValue("delta");
      model_pointer->set_model_parameter("delta", delta);
    }
    if (not tau_fixed and not DT_fixed) {
      double tau = getParameterValue("tau");
      model_pointer->set_model_parameter("tau", tau);
      if (tau < 1e-10)
        no_T = true;
    }
    if (not lambda_fixed) {
      double lambda = getParameterValue("lambda");
      model_pointer->set_model_parameter("lambda", lambda);
    }
    if (DT_fixed) {
      double tau = getParameterValue("tau");
      model_pointer->set_model_parameter("tau", tau);
      double delta = tau * model_pointer->scalar_parameter["DT_ratio"];
      model_pointer->set_model_parameter("delta", delta);
    }
    if (MLOR) {
      double O_R = getParameterValue("O_R");
      model_pointer->set_model_parameter("O_R", O_R);
    }

    for (int i = 0; i < public_ml_branch_ids.size(); i++) {
      int e = public_ml_branch_ids[i];
      string name = public_ml_ratetype_names[i];
      stringstream branch;
      branch << e;
      scalar_type multiplier =
          getParameterValue("rm_" + name + "_" + branch.str());
      model_pointer->vector_parameter[name][e] = multiplier;
    }

    model_pointer->calculate_undatedEs();
    double y = -log(model_pointer->pun(ale_pointer, false, no_T));
    // cout <<endl<< "delta=" << model_pointer->scalar_parameter["delta_avg"] <<
    // "\t tau=" <<model_pointer->scalar_parameter["tau_avg"] << "\t lambda=" <<
    // model_pointer->scalar_parameter["lambda_avg"] << "\t ll="
    //           << -y <<endl;
    fval_ = y;
  }
};

int main(int argc, char **argv) {
  cout << "ALEml_undated using ALE v" << ALE_VERSION << endl;
  bool no_T = false;
  if (argc < 3) {
    cout << "\nUsage:\n ./ALEml_undated species_tree.newick "
            "gene_tree_sample.ale sample=number_of_samples seed=integer "
            "separators=gene_name_separator O_R=OriginationAtRoot "
            "delta=DuplicationRate tau=TransferRate lambda=LossRate "
            "beta=weight_of_sequence_evidence "
            "fraction_missing=file_with_fraction_of_missing_genes_per_species "
            "output_species_tree=n S_branch_lengths:root_length "
            "rate_multiplier:rate_name:branch_id:value"
         << endl;
    cout << "\n1st example: we fix the DTL values and do not perform any "
            "optimization \n ./ALEml_undated species_tree.newick "
            "gene_tree_sample.ale sample=100 separators=_ delta=0.05 tau=0.1 "
            "lambda=0.2 "
         << endl;
    cout << "\n2nd example: we fix the T value to 0 to get a DL-only model and "
            "optimize the DL parameters \n ./ALEml_undated species_tree.newick "
            "gene_tree_sample.ale sample=100 separators=_ tau=0\n"
         << endl;
    cout << "\n3rd example: we fix the ratio of D and T rates to D/T=1/20\n "
            "./ALEml_undated species_tree.newick gene_tree_sample.ale "
            "sample=100 separators=_ DT=0.05\n"
         << endl;
    cout << "\n4th example: we provide a file giving the expected fraction of "
            "missing genes in each species \n ./ALEml_undated "
            "species_tree.newick gene_tree_sample.ale sample=100 separators=_ "
            "fraction_missing=fraction_missing.txt\n"
         << endl;
    cout << "\n5th example: same as 3rd, but outputs the annotated species "
            "tree to a file \n ./ALEml_undated species_tree.newick "
            "gene_tree_sample.ale sample=100 separators=_ "
            "fraction_missing=fraction_missing.txt output_species_tree=y\n"
         << endl;
    cout << "\n6th example: use species tree branch lengths as fixed rate "
            "multipliers with root length specifed as 0.2 (default 1.)\n "
            "./ALEml_undated species_tree.newick gene_tree_sample.ale "
            "S_branch_lengths:0.2 \n"
         << endl;
    cout << "\n6.1th example: use fixed branchrate multiplier for rate of Ts "
            "_to_ branch 43 with value 0.0 (i.e. no transfer to branch)\n "
            "./ALEml_undated species_tree.newick gene_tree_sample.ale "
            "rate_multiplier:tau_to:43:0.0 \n"
         << endl;
    cout << "\n6.2th example: use fixed branchrate multiplier for rate of Ts "
            "_from_ branch 43 with value 0.0 (i.e no transfer from branch)\n "
            "./ALEml_undated species_tree.newick gene_tree_sample.ale "
            "rate_multiplier:tau_from:43:0.0 \n"
         << endl;
    cout << "\n6.3th example: use fixed branchrate multiplier for rate Ds on "
            "branch 43 with value 0.0 (no duplications on branch)\n "
            "./ALEml_undated species_tree.newick gene_tree_sample.ale "
            "rate_≈repl:delta:43:0.0 \n"
         << endl;
    cout << "\n6.4th example: use fixed branchrate multiplier for rate Ls on "
            "branch 43 with value 0.0 (no losses on branch)\n ./ALEml_undated "
            "species_tree.newick gene_tree_sample.ale "
            "rate_multiplier:lambda:43:0.0 \n"
         << endl;
    cout << "\n6.1.1 example: optimize branchrate multiplier for rate of Ts "
            "_to_ branch 43 (same syntax for all other multipliers!)\n "
            "./ALEml_undated species_tree.newick gene_tree_sample.ale "
            "S_branch_lengths:-2 \n"
         << endl;
    cout << "\n6.5th example: use fixed multiplier for Origination on branch "
            "43 with value 0.0 (no origination on branch)\n ./ALEml_undated "
            "species_tree.newick gene_tree_sample.ale rate_multiplier:O:43:0.0 "
            "\n"
         << endl;
    cout
        << "\n6.6th example: fixed multiplier for Origination on branch 43 to "
           "1 and zero elsewhere (no origination outside 43)\n ./ALEml_undated "
           "species_tree.newick gene_tree_sample.ale rate_multiplier:O:43:-1 \n"
        << endl;

    return 0;
  }

  // we need a rooted species tree in newick format
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
  bool DT_fixed = false;

  scalar_type delta = 1e-2, tau = 1e-2, lambda = 1e-1, DT_ratio = 0.05;
  string fractionMissingFile = "";
  string outputSpeciesTree = "";
  map<string, map<int, scalar_type>> rate_multipliers;
  vector<int> ml_branch_ids;
  vector<string> ml_ratetype_names;
  bool MRP = false;
  bool MLOR = false;

  model->set_model_parameter("undatedBL", false);
  model->set_model_parameter("reldate", false);
  MLOR = false;

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
      if (tau < 1e-10) {
        no_T = true;
        cout << "# tau fixed to no transfer!" << endl;
        tau = 1e-19;
      } else
        cout << "# tau fixed to " << tau << endl;
    } else if (tokens[0] == "lambda") {
      lambda = atof(tokens[1].c_str());
      lambda_fixed = true;
      cout << "# lambda fixed to " << lambda << endl;
    } else if (tokens[0] == "DT") {
      DT_ratio = atof(tokens[1].c_str());
      DT_fixed = true;
      model->set_model_parameter("DT_ratio", DT_ratio);
      cout << "# D/T ratio fixed to " << model->scalar_parameter["DT_ratio"]
           << endl;
    } else if (tokens[0] == "O_R") {
      O_R = atof(tokens[1].c_str());
      cout << "# O_R set to " << O_R << endl;
    } else if (tokens[0] == "beta") {
      beta = atof(tokens[1].c_str());
      cout << "# beta set to " << beta << endl;
    } else if (tokens[0] == "fraction_missing") {
      fractionMissingFile = tokens[1];
      cout << "# File containing fractions of missing genes set to "
           << fractionMissingFile << endl;
    } else if (tokens[0] == "S_branch_lengths") {
      model->set_model_parameter("undatedBL", true);
      if (tokens.size() == 1) {
        model->set_model_parameter("root_BL", 1);
        cout << "# unsing branch lengths of input S tree as rate multipliers "
                "with 1 at root! "
             << endl;
      } else {
        scalar_type root_rm = atof(tokens[1].c_str());
        model->set_model_parameter("root_BL", root_rm);
        cout << "# unsing branch lengths of input S tree as rate multipliers "
                "with "
             << root_rm << " at root! " << endl;
      }

    } else if (tokens[0] == "reldate") {
      cout << "Respecting realtive ages from input S tree, please make sure "
              "input S tree is ultrametric!"
           << endl;
      model->set_model_parameter("reldate", true);
    } else if (tokens[0] == "MLOR") {
      cout << "Optimizing root origination multiplier." << endl;
      MLOR = true;
    }

    else if (tokens[0] == "rate_multiplier") {
      string rate_name = tokens[1];
      int e = atoi(tokens[2].c_str());
      scalar_type rm = atof(tokens[3].c_str());
      if (rm >= -1) {
        cout << "# rate multiplier for rate " << rate_name
             << " on branch with ID " << e << " set to " << rm << endl;
        rate_multipliers["rate_multiplier_" + rate_name][e] = rm;
      } else {
        cout << "# rate multiplier for rate " << rate_name
             << " on branch with ID " << e << " to be optimized " << endl;
        ml_branch_ids.push_back(e);
        ml_ratetype_names.push_back("rate_multiplier_" + rate_name);
      }
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

    } else if (tokens[0] == "seed") {
      long seed = atoi(tokens[1].c_str());
      cout << "Set random seed to " << seed << endl;
      RandomTools::setSeed(seed);
    } else if (tokens[0] == "MRP") {
      bool MRP = true;
    }
  }

  model->set_model_parameter("BOOTSTRAP_LABELS", "yes");
  model->construct_undated(Sstring, fractionMissingFile);

  for (auto it = rate_multipliers.begin(); it != rate_multipliers.end(); it++)
    for (auto jt = (*it).second.begin(); jt != (*it).second.end(); jt++) {
      model->vector_parameter[(*it).first][(*jt).first] = (*jt).second;
    }

  model->set_model_parameter("seq_beta", beta);
  model->set_model_parameter("O_R", O_R);
  // a set of inital rates
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);

  // calculate_EGb() must always be called after changing rates to calculate E-s
  // and G-s cf. http://arxiv.org/abs/1211.4606
  model->calculate_undatedEs();
  cout << "Reconciliation model initialised, starting DTL rate optimisation"
       << ".." << endl;

  // we use the Nelder–Mead method implemented in Bio++

  Function *f =
      new p_fun(model, ale, ml_branch_ids, ml_ratetype_names, delta, tau,
                lambda, delta_fixed, tau_fixed, lambda_fixed, DT_fixed, MLOR);

  // ReparametrizationFunctionWrapper rpf(f, false);
  // ThreePointsNumericalDerivative tpnd(&rpf);
  // tpnd.setParametersToDerivate(rpf.getParameters().getParameterNames());

  Optimizer *optimizer;
  if (delta_fixed or tau_fixed or lambda_fixed or DT_fixed) {
    optimizer = new DownhillSimplexMethod(f);
  } else {
    optimizer = new DownhillSimplexMethod(f);
  }
  optimizer->setProfiler(0);
  optimizer->setMessageHandler(0);
  optimizer->setVerbose(2);

  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  optimizer->init(f->getParameters()); // Here we optimize all parameters, and
                                       // start with the default values.

  scalar_type mlll;

  if (not(delta_fixed and tau_fixed and lambda_fixed and
          not MLOR)) // not all rates fixed
  {
    cout << "#optimizing rates" << endl;
    optimizer->optimize();
    if (not delta_fixed and not DT_fixed)
      delta = optimizer->getParameterValue("delta");
    if (not tau_fixed)
      tau = optimizer->getParameterValue("tau");
    if (DT_fixed)
      delta = tau * model->scalar_parameter["DT_ratio"];
    if (not lambda_fixed)
      lambda = optimizer->getParameterValue("lambda");
    if (MLOR)
      O_R = optimizer->getParameterValue("O_R");
    mlll = -optimizer->getFunctionValue();

  } else {
    mlll = log(model->pun(ale, false, no_T));
  }
  stringstream ml_rate_multipliers;
  for (int i = 0; i < ml_branch_ids.size(); i++) {
    int e = ml_branch_ids[i];
    string name = ml_ratetype_names[i];
    stringstream branch;
    branch << e;
    scalar_type multiplier =
        optimizer->getParameterValue("rm_" + name + "_" + branch.str());
    ml_rate_multipliers << name << "\t" << e << "\t" << multiplier << ";\n";
  }

  cout << endl
       << "ML rates: "
       << " delta=" << delta << "; tau=" << tau << "; lambda=" << lambda
       << "; O_R=" << O_R //<<"; sigma="<<sigma_hat
       << "." << endl;

  if (ml_branch_ids.size() > 0)
    cout << "ML rate multipliers:\n" << ml_rate_multipliers.str();

  cout << "LL=" << mlll << endl;

  cout << "Sampling reconciled gene trees.." << endl;
  vector<string> sample_strings;
  vector<Tree *> sample_trees;
  boost::progress_display pd(samples);

  for (int i = 0; i < samples; i++) {
    ++pd;
    string sample_tree = model->sample_undated(no_T);
    sample_strings.push_back(sample_tree);
    if (ale->last_leafset_id > 3 and MRP) {

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
  string outname = ale_file + ".uml_rec";
  ofstream fout(outname.c_str());
  fout << "#ALEml_undated using ALE v" << ALE_VERSION
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
  if (ml_branch_ids.size() > 0) {
    fout << "ML rate multipliers:\n" << ml_rate_multipliers.str();
    fout << endl;
  }
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
  fout << "# of\t "
          "Duplications\tTransfers\tLosses\tOriginations\tcopies\tsingletons\te"
          "xtinction_prob\tpresence\tLL"
       << endl;
  fout << model->counts_string_undated(samples);
  fout.close();

  // Outputting the species tree to its own file:
  if (outputSpeciesTree != "") {
    ofstream fout(outputSpeciesTree.c_str());
    fout << model->string_parameter["S_with_ranks"] << endl;
    fout.close();
  }

  cout << "Results in: " << outname << endl;
  if (ale->last_leafset_id > 3 and MRP) {
    //      cout << "Calculating consensus tree."<<endl;
    // cout << TreeTools::treeToParenthesis(sample_trees[0])<<std::endl;

    cout << "Calculating MRP tree." << endl;
    Tree *con_tree = TreeTools::MRP(sample_trees);
    // cout << "Calculating threshold consensus tree."<<endl;

    // con_tree= TreeTools::thresholdConsensus(sample_trees,0.5);

    string con_name = ale_file + ".ucons_tree";

    ofstream con_out(con_name.c_str());
    con_out << "#ALEsample using ALE v" << ALE_VERSION
            << " by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;" << endl;
    // TreeTools::computeBootstrapValues(*con_tree,sample_trees);
    string con_tree_sup = TreeTemplateTools::treeToParenthesis(*con_tree);
    con_out << con_tree_sup << endl;
    cout << endl << "Consensus tree in " << con_name << endl;
  }

  string t_name = ale_file + ".uTs";
  ofstream tout(t_name.c_str());
  tout << "#from\tto\tfreq.\n";

  for (int e = 0; e < model->last_branch; e++)
    for (int f = 0; f < model->last_branch; f++)
      if (model->T_to_from[e][f] > 0) {
        if (e < model->last_leaf)
          tout << model->node_name[model->id_nodes[e]] << "(" << e << ")";
        else
          tout << "\t" << e;
        if (f < model->last_leaf)
          tout << "\t" << model->node_name[model->id_nodes[f]] << "(" << f
               << ")";
        else
          tout << "\t" << f;
        tout << "\t" << model->T_to_from[e][f] / samples << endl;
      }
  tout.close();
  cout << "Transfers in: " << t_name << endl;
  return 0;
}
