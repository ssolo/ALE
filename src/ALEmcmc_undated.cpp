
#include "exODT.h"
#include "ALE_util.h"

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>

using namespace std;
using namespace bpp;


/************************************************************************
 * Scaling move to change the value of a Real Positive (e.g. a rate).
 * The lambda parameter here represents a scale.
 ************************************************************************/
//This version of the scaling function ensures that the new value is not larger than some maximum value
double scaleDoubleConstrained ( const double& value, const double& maxi, const double& lambda, double& hastingsRatio, const bool verbose ) {
	double newValue = value;
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = RandomTools::giveRandomNumberBetweenZeroAndEntry ( 1.0 );
    double scalingFactor = std::exp( lambda * ( u - 0.5 ) );
    newValue *= scalingFactor;
    if (newValue < 0.00001 ) {
		newValue = 0.00001;
	}
	if (newValue > maxi ) {
		newValue = maxi;
	}

    // compute the Hastings ratio
    hastingsRatio = scalingFactor ;
	return newValue;
}

/************************************************************************
 * Exponential probability for a double (e.g. a rate),
 * given the parameter of the exponential distribution (a rate = 1/scale).
 * P(value) = lambda exp (-lambda*value)
 * mean:1/lambda
 ************************************************************************/
double computeExponentialLogProbability ( const double& param, const double& value ) {
	return (std::log(param) - param * value);
}


//////////////////////////////////////////
//////////////////////////////////////////

double computeLogLk (exODT_model* model, approx_posterior * ale, const double &o, const double &d, const double &t, const double &l) {
  model->set_model_parameter("O_R", o);
  model->set_model_parameter("delta", d);
  model->set_model_parameter("tau", t);
  model->set_model_parameter("lambda", l);
  //calculate_EGb() must always be called after changing rates to calculate E-s and G-s
  //cf. http://arxiv.org/abs/1211.4606
  model->calculate_undatedEs();
  double ll = log(model->pun(ale));
  return ll;
}

//////////////////////////////////////////
//////////////////////////////////////////

double computeLogPrior ( const double &o, const double &d, const double &t, const double &l, const double &priorOrigination, const double &priorDelta, const double &priorTau,  const double &priorLambda) {
  double pp = 0.0;
  pp += computeExponentialLogProbability ( priorOrigination, o);
  pp += computeExponentialLogProbability ( priorDelta, d);
  pp += computeExponentialLogProbability ( priorTau, t);
  pp += computeExponentialLogProbability ( priorLambda, l);
  return pp;
}

//////////////////////////////////////////
//////////////////////////////////////////

void     acceptMove(double &currentOrigination, double &currentDelta, double &currentTau, double &currentLambda, const double &newOrigination, const double &newDelta, const double &newTau, const double &newLambda) {
  currentOrigination = newOrigination;
  currentDelta = newDelta;
  currentTau = newTau;
  currentLambda = newLambda;
}

//////////////////////////////////////////
//////////////////////////////////////////

void     rejectMove( const double &currentOrigination, const double &currentDelta, const double &currentTau, const double &currentLambda, double &newOrigination, double &newDelta, double &newTau, double &newLambda) {
newOrigination = currentOrigination;
newDelta = currentDelta;
newTau = currentTau;
newLambda = currentLambda;
}

//////////////////////////////////////////
//////////////////////////////////////////

void sampleTree (exODT_model* model, approx_posterior * ale, double &o, double &d, double &t, double &l, vector <string> &sample_strings, vector <Tree*> &sample_trees) {
  model->set_model_parameter("O_R", o);
  model->set_model_parameter("delta", d);
  model->set_model_parameter("tau", t);
  model->set_model_parameter("lambda", l);
  //calculate_EGb() must always be called after changing rates to calculate E-s and G-s
  //cf. http://arxiv.org/abs/1211.4606
  model->calculate_undatedEs();
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


void fillTToFrom(exODT_model* model, map<string, double>& tToFrom) {
    for (int e=0;e<model->last_branch;e++) {
      for (int f=0;f<model->last_branch;f++) {
        if  (model->T_to_from[e][f]>0)
  	{
      string name1, name2, names;
  	  if (e<model->last_leaf) {
  	    name1 = model->node_name[model->id_nodes[e]];
      }
  	  else {
  	    name1 = std::to_string(e);
      }
  	  if (f<model->last_leaf) {
  	    name2 = model->node_name[model->id_nodes[f]];
      }
  	  else {
  	    name2 = std::to_string(f);
      }
      names = name1 + "\t" + name2;
      map<string, double>::iterator it = tToFrom.find(names);
      if(it != tToFrom.end())
      {
        //element found;
        tToFrom[names] += model->T_to_from[e][f];
      }
      else {
        tToFrom[names] = model->T_to_from[e][f];
      }
  	}
  }
  }
  return;
}


//////////////////////////////////////////
//////////////////////////////////////////

int main(int argc, char ** argv)
{
  cout << "ALEmcmc using ALE v"<< ALE_VERSION <<endl;

  if (argc<3)
  {
    cout << "\nUsage:\n ./ALEmcmc_undated species_tree.newick gene_tree_sample.ale sample=number_of_samples separators=gene_name_separator O_R=OriginationAtRootPrior delta=DuplicationRatePrior tau=TransferRatePrior lambda=LossRatePrior sampling_rate=sampling_rate beta=weight_of_sequence_evidence fraction_missing=file_with_fraction_of_missing_genes_per_species output_species_tree=n" << endl;
    cout << "\nExample:\n ./ALEmcmc_undated species_tree.newick gene_tree_sample.ale sample=100 separators=_ O_R=1 delta=0.01 tau=0.01 lambda=0.1 sampling_rate=10 beta=1\n" << endl;
		cout << "\n2nd example: we provide a file giving the expected fraction of missing genes in each species \n ./ALEmcmc_undated species_tree.newick gene_tree_sample.ale sample=100 separators=_ fraction_missing=fraction_missing.txt\n" << endl;
		cout << "\n3rd example: same as 2nd, but outputs the annotated species tree to a file \n ./ALEmcmc_undated species_tree.newick gene_tree_sample.ale sample=100 separators=_ fraction_missing=fraction_missing.txt output_species_tree=y\n" << endl;

    return 0;
  }

  //we need a rooted species tree in newick format
  string Sstring;
  string S_treefile=argv[1];
	if (!fexists(argv[1])) {
		cout << "Error, file "<<argv[1] << " does not seem accessible." << endl;
		exit(1);
	}
  ifstream file_stream_S (argv[1]);
  getline (file_stream_S,Sstring);
  cout << "Read species tree from: " << argv[1] <<".."<<endl;
  //we need an .ale file containing observed conditional clade probabilities
  //cf. ALEobserve
  string ale_file=argv[2];
  approx_posterior * ale;
  ale=load_ALE_from_file(ale_file);
  cout << "Read summary of tree sample for "<<ale->observations<<" trees from: " << ale_file <<".."<<endl;

	//Getting the radical for output files:
  vector<string> tokens;
  boost::split(tokens,ale_file,boost::is_any_of("/"),boost::token_compress_on);
  ale_file=tokens[tokens.size()-1];


  //we initialise a coarse grained reconciliation model for calculating the sum
  exODT_model* model=new exODT_model();

  scalar_type samples=100;

	//a set of inital rates
  double priorOrigination =1.0,  priorDelta=0.01,priorTau=0.01,priorLambda=0.1;
	size_t sampling_rate = 10;
	scalar_type beta=1;

	string fractionMissingFile = "";
	string outputSpeciesTree = "";

	for (int i=3;i<argc;i++)
	{
		string next_field=argv[i];
		vector <string> tokens;
		boost::split(tokens,next_field,boost::is_any_of("="),boost::token_compress_on);
		if (tokens[0]=="sample")
		samples=atof(tokens[1].c_str());
		else if (tokens[0]=="separators")
		model->set_model_parameter("gene_name_separators", tokens[1]);
		else if (tokens[0]=="delta")
		{
			priorDelta=atof(tokens[1].c_str());
			cout << "# priorDelta fixed to " << priorDelta << endl;
		}
		else if (tokens[0]=="tau")
		{
			priorTau=atof(tokens[1].c_str());
			cout << "# priorTau fixed to " << priorTau << endl;
		}
		else if (tokens[0]=="lambda")
		{
			priorLambda=atof(tokens[1].c_str());
			priorLambda=true;
			cout << "# priorLambda fixed to " << priorLambda << endl;

		}
		else if (tokens[0]=="O_R")
		{
			priorOrigination=atof(tokens[1].c_str());
			cout << "# priorOrigination set to " << priorOrigination << endl;
		}
		else if (tokens[0]=="beta")
		{
			beta=atof(tokens[1].c_str());
			cout << "# beta set to " << beta << endl;
		}
		else if (tokens[0]=="sampling_rate")
		{
			sampling_rate=atoi(tokens[1].c_str());
			cout << "# sampling_rate set to " << sampling_rate << endl;
		}
		else if (tokens[0]=="fraction_missing")
		{
			fractionMissingFile=tokens[1];
			cout << "# File containing fractions of missing genes set to " << fractionMissingFile << endl;
		}
		else if (tokens[0]=="output_species_tree")
		{
			std::string valArg = boost::algorithm::to_lower_copy(tokens[1]);
			if (valArg == "y" || valArg == "ye" || valArg == "yes" ) {
					outputSpeciesTree= ale_file + ".spTree";
					cout << "# outputting the annotated species tree to "<< outputSpeciesTree << endl;
			}
		}

	}


  model->set_model_parameter("BOOTSTRAP_LABELS","yes");
	model->set_model_parameter("seq_beta", beta);

  model->construct_undated(Sstring, fractionMissingFile);

  double currentOrigination = RandomTools::randExponential(priorOrigination) ;
  double currentDelta = RandomTools::randExponential(priorDelta) ;
  double currentTau = RandomTools::randExponential(priorTau) ;
  double currentLambda = RandomTools::randExponential(priorLambda) ;

  double newOrigination = currentOrigination;
  double newDelta = currentDelta;
  double newTau = currentTau;
  double newLambda = currentLambda;

  double currentLogLikelihood = computeLogLk (model, ale, currentOrigination, currentDelta, currentTau, currentLambda);
  double newLogLikelihood = currentLogLikelihood;

  double currentLogPrior = computeLogPrior (currentOrigination, currentDelta, currentTau, currentLambda, priorOrigination, priorDelta, priorTau, priorLambda) ;
  double newLogPrior = currentLogPrior;

  std::cout << "Initial logLK: "<< currentLogLikelihood << " and logPrior: "<< currentLogPrior <<std::endl;

  cout << "Reconciliation model initialised, starting DTL rate sampling" <<".."<<endl;


  size_t originationId = 0;
  size_t deltaId = 1;
  size_t lambdaId = 2;
  size_t tauId = 3;

  std::vector<double> moveWeights;
  moveWeights.push_back( 1 ) ; // originationId
  moveWeights.push_back( 1 ); // deltaId
  moveWeights.push_back( 1 ); // lambdaId
  moveWeights.push_back( 1 ); // tauId

  size_t move = 0;
  std::vector< int > order = VectorTools::seq ( 0, (int) moveWeights.size(), 1 );
  double maxSumDTL = 10;
  double maxOrigination = 1000000;
  double hastingsRatio = 0.0;
  bool verbose = false;
  double scale = 1;
  std::vector<double> scaleMoveParameters ;
  scaleMoveParameters.push_back (0.1) ;
  scaleMoveParameters.push_back (1) ;
  scaleMoveParameters.push_back (10) ;
  std::vector<double> scaleWeights;
  scaleWeights.push_back(1);
  scaleWeights.push_back(1);
  scaleWeights.push_back(1);
  double threshold = 0.0;
  double acceptanceProbability = 0.0;
  std::vector<int> orderScaleMoveParameters = VectorTools::seq ( 0, (int) scaleMoveParameters.size(), 1 );
  vector <string> sample_strings;
  vector <Tree*> sample_trees;
  size_t i = 0;
  // Summary variables
  double numSpeciations = 0.0;
  double numDuplications = 0.0;
  double numTransfers = 0.0;
  double numLosses = 0.0;
  map<string, double> tToFrom ;
  //boost::progress_display pd( i );

  string mcmcoutname=ale_file+"_umcmc.csv";
  ofstream mcmcout( mcmcoutname.c_str() );

  mcmcout << "Iteration" << "\t"<< "LogLk" << "\t" << "LogPrior"  << "\t" << "Origination" << "\t" << "Delta" << "\t" << "Tau" << "\t" << "Lambda" <<std::endl;
  std::cout << "Iteration" << "\t"<< "LogLk" << "\t" << "LogPrior"  << "\t" << "Origination" << "\t" << "Delta" << "\t" << "Tau" << "\t" << "Lambda" <<std::endl;

  // BURNIN loop
  size_t burninLength = 100;
  std::cout << "BURNIN during "<<burninLength<<" iterations."<<std::endl;
  std::cout << "LogLk" << "\t" << "LogPrior" << "\t" << "Origination" << "\t" << "Delta" << "\t" << "Tau" << "\t" << "Lambda"<<std::endl;
  for (i = 0 ; i < burninLength ; ++i) {
    move = RandomTools::pickOne( order, moveWeights, true );
    scale = scaleMoveParameters[RandomTools::pickOne( orderScaleMoveParameters, scaleWeights, true )];
    if(move == originationId) {
      newOrigination = scaleDoubleConstrained ( currentOrigination, maxOrigination, scale, hastingsRatio, verbose ) ;
    }
    else if (move == deltaId) {
      newDelta = scaleDoubleConstrained ( currentDelta, maxSumDTL - currentLambda - currentTau, scale, hastingsRatio, verbose ) ;
    }
    else if (move == lambdaId) {
      newLambda = scaleDoubleConstrained ( currentLambda, maxSumDTL - currentDelta - currentTau, scale, hastingsRatio, verbose ) ;
    }
    else if (move == tauId) {
      newTau = scaleDoubleConstrained ( currentTau, maxSumDTL - currentLambda - currentDelta, scale, hastingsRatio, verbose ) ;
    }
    newLogLikelihood = computeLogLk(model, ale, newOrigination, newDelta, newTau, newLambda);
    newLogPrior = computeLogPrior(newOrigination, newDelta, newTau, newLambda, priorOrigination, priorDelta, priorTau, priorLambda);
    //Accept or reject?
    acceptanceProbability = exp (  ( newLogLikelihood + newLogPrior ) -  ( currentLogLikelihood + currentLogPrior) ) * hastingsRatio ;
    threshold = RandomTools::giveRandomNumberBetweenZeroAndEntry ( 1.0 );
    if (acceptanceProbability > threshold) { //accept
      acceptMove(currentOrigination, currentDelta, currentTau, currentLambda, newOrigination, newDelta, newTau, newLambda);
      currentLogLikelihood = newLogLikelihood;
      currentLogPrior = newLogPrior;
    }
    else {
      rejectMove(currentOrigination, currentDelta, currentTau, currentLambda, newOrigination, newDelta, newTau, newLambda);
    }
    std::cout << i <<"\t"<< currentLogLikelihood << "\t" << currentLogPrior  << "\t" << currentOrigination << "\t" << currentDelta << "\t" << currentTau << "\t" << currentLambda <<std::endl;
  }


  // MCMC loop
  std::cout << "MCMC during "<<samples * sampling_rate<<" iterations."<<std::endl;
  std::cout << "LogLk" << "\t" << "LogPrior" << "\t" << "Origination" << "\t" << "Delta" << "\t" << "Tau" << "\t" << "Lambda"<<std::endl;
  for (i = 0 ; i < samples * sampling_rate ; ++i) {
    move = RandomTools::pickOne( order, moveWeights, true );
    scale = scaleMoveParameters[RandomTools::pickOne( orderScaleMoveParameters, scaleWeights, true )];
    if(move == originationId) {
      newOrigination = scaleDoubleConstrained ( currentOrigination, maxOrigination, scale, hastingsRatio, verbose ) ;
    }
    else if (move == deltaId) {
      newDelta = scaleDoubleConstrained ( currentDelta, maxSumDTL - currentLambda - currentTau, scale, hastingsRatio, verbose ) ;
    }
    else if (move == lambdaId) {
      newLambda = scaleDoubleConstrained ( currentLambda, maxSumDTL - currentDelta - currentTau, scale, hastingsRatio, verbose ) ;
    }
    else if (move == tauId) {
      newTau = scaleDoubleConstrained ( currentTau, maxSumDTL - currentLambda - currentDelta, scale, hastingsRatio, verbose ) ;
    }
    newLogLikelihood = computeLogLk(model, ale, newOrigination, newDelta, newTau, newLambda);
    newLogPrior = computeLogPrior(newOrigination, newDelta, newTau, newLambda, priorOrigination, priorDelta, priorTau, priorLambda);
    //Accept or reject?
    acceptanceProbability = exp (  ( newLogLikelihood + newLogPrior ) -  ( currentLogLikelihood + currentLogPrior) ) * hastingsRatio ;
    threshold = RandomTools::giveRandomNumberBetweenZeroAndEntry ( 1.0 );
    if (acceptanceProbability > threshold) { //accept
      acceptMove(currentOrigination, currentDelta, currentTau, currentLambda, newOrigination, newDelta, newTau, newLambda);
      currentLogLikelihood = newLogLikelihood;
      currentLogPrior = newLogPrior;
    }
    else {
      rejectMove(currentOrigination, currentDelta, currentTau, currentLambda, newOrigination, newDelta, newTau, newLambda);
    }
    if (i % sampling_rate == 0 ) {
      model->MLRec_events.clear();
      model->reset_T_to_from();
      sampleTree (model, ale, currentOrigination, currentDelta, currentTau, currentLambda, sample_strings, sample_trees);
      numSpeciations += model->MLRec_events["S"];
      numDuplications += model->MLRec_events["D"];
      numTransfers += model->MLRec_events["T"];
      numLosses += model->MLRec_events["L"];
      fillTToFrom(model, tToFrom);
      std::cout << i <<"\t"<< currentLogLikelihood << "\t" << currentLogPrior  << "\t" << currentOrigination << "\t" << currentDelta << "\t" << currentTau << "\t" << currentLambda <<std::endl;
    }
    mcmcout << i <<"\t"<< currentLogLikelihood << "\t" << currentLogPrior  << "\t" << currentOrigination << "\t" << currentDelta << "\t" << currentTau << "\t" << currentLambda <<std::endl;
  }
  mcmcout.close();

  string outname=ale_file+".umcmc_rec";
  ofstream fout( outname.c_str() );
  fout <<  "#ALEmcmc_undated using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl<<endl;
  fout << "S:\t"<<model->string_parameter["S_with_ranks"] <<endl;
  fout << endl;
  fout << "Input ale from:\t"<<ale_file<<endl;
  //fout << ">logl: " << mlll << endl;
  //fout << "rate of\t Duplications\tTransfers\tLosses" <<endl;
  //fout << "ML \t"<< delta << "\t" << tau << "\t" << lambda //<< "'t" << sigma_hat << endl;
  fout << endl;
  fout << samples << " reconciled G-s:\n"<<endl;
  for (int i=0;i<samples;i++)
    {
      fout<<sample_strings[i]<<endl;
    }

  fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" <<endl;
  fout <<"Total \t"<< numDuplications/samples << "\t" << numTransfers/samples << "\t" << numLosses/samples<< "\t" << numSpeciations/samples <<endl;
  fout << endl;
  fout.close();
  /*fout << "# of\t Duplications\tTransfers\tLosses\tcopies" <<endl;
  fout << model->counts_string_undated(samples);
*/

// Outputting the species tree to its own file:
  if (outputSpeciesTree != "") {
    ofstream fout( outputSpeciesTree.c_str() );
    fout <<model->string_parameter["S_with_ranks"] <<endl;
    fout.close();
  }

  cout << "Results in: " << outname << endl;
  if (ale->last_leafset_id>3)
    {
      cout << "Calculating consensus tree."<<endl;
      Tree* con_tree= TreeTools::thresholdConsensus(sample_trees,0.5);

      string con_name=ale_file+".ucons_mcmc_tree";

      ofstream con_out( con_name.c_str() );
      con_out <<  "#ALEmcmc_undated using ALE v"<< ALE_VERSION <<" by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"<<endl;
      TreeTools::computeBootstrapValues(*con_tree,sample_trees);
      string con_tree_sup=TreeTemplateTools::treeToParenthesis(*con_tree);
      con_out << con_tree_sup << endl;
      cout << endl<< "Consensus tree in " << con_name<< endl;
			con_out.close();
    }

  string t_name=ale_file+"_mcmc.uTs";
  ofstream tout( t_name.c_str() );
  tout <<"#from\tto\tfreq.\n";
  std::map<std::string, double>::const_iterator it;
  for(it = tToFrom.begin(); it != tToFrom.end(); it++) {
     tout << "\t" << it->first << "\t" <<  it->second/samples << endl;
  }
	tout.close();
  cout << "Transfers in: " << t_name << endl;
  return 0;
}
