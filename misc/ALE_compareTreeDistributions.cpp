#include "ALE.h"
#include "ALE_util.h"
using namespace std;

//Compilation:  g++ -o ALE_compareTreeDistributions ALE_compareTreeDistributions.cpp ALE.h ALE.cpp ALE_util.h ALE_util.cpp exODT.h exODT.cpp model_omp.cpp -std=c++0x -I/usr/local/include  -L. -L/usr/local/lib  -lbpp-core -lbpp-seq -lbpp-phyl
//  g++ -g -o ALE_compareTreeDistributions ALE_compareTreeDistributions.cpp ALE.h ALE.cpp ALE_util.h ALE_util.cpp exODT.h exODT.cpp model_omp.cpp -std=c++0x -I/usr/include  -I/home/ssolo/newest_bpp/include -L. -L/home/ssolo/newest_bpp/lib -L/usr/lib  -lbpp-core -lbpp-seq -lbpp-phyl

int main(int argc, char ** argv)
{
  if (argc == 1 ) {
	std::cout << "\tUsage: ALE_compareTreeDistributions treeDist1.trees burnin1 treeDist2.trees burnin2\n" <<std::endl;
  }
else {
  //First tree distribution
  string ale_file=argv[1];
  string ale_name=ale_file+".ale";
  approx_posterior * ale;
  vector<int> burnin;

  if (argc>2) 
    burnin.push_back(atoi(argv[2]));  
	vector<int> every ;
if (argc>3) 
    every.push_back(atoi(argv[3]));

  ale=observe_ALE_from_file(ale_file, burnin[0], every[0]);

  cout << "# observe "<< ale->observations << "trees from: " <<  argv[1] << endl;
  ale->save_state(ale_name);
  cout << "# saved in "<< ale_name<<endl;

  //Second tree distribution
  string ale_file2=argv[4];
  string ale_name2=ale_file+".ale";
  approx_posterior * ale2;

  if (argc>4) 
    burnin.push_back(atoi(argv[5]));
	int every2 = 0;
  if (argc>5) 
    every.push_back(atoi(argv[6]));

  ale2=observe_ALE_from_file(ale_file2, burnin[1], every[1]);
  cout << "# observe "<< ale2->observations << "trees from: " <<  argv[4] << endl;
  ale2->save_state(ale_name2);
  cout << "# saved in "<< ale_name2<<endl;

  //Now we want to compute the probabilities of all trees included in the files, according to the two ales.  
  //First, we get all trees and put them in a single vector

  vector<string> trees;
  std::vector< string > fnames;
  fnames.push_back(argv[1]);
  fnames.push_back(argv[4]);
  size_t fileId = 0;
  for (vector<string>::iterator it=fnames.begin();it!=fnames.end();it++)
    {
      string fname=(*it);
      ifstream file_stream (fname.c_str());
      int tree_i=0;  
      if (file_stream.is_open())  //  ########## read trees ############
	{
	  while (! file_stream.eof())
	    {
	      string line;
	      getline (file_stream,line);
	      if (line.find("(")!=line.npos )
		{
		  tree_i++;		 
		  if (tree_i%every[fileId]==0 ) trees.push_back(line);			     
		}
	    }
	}
    fileId++;
    }

  string outFile;
   if (argc>4) 
    outFile=atoi(argv[5]);
  else
	outFile = "outFile.txt";

  ofstream myfile;
  myfile.open (outFile);
  myfile << "ALE1_proba\tALE2_proba"<<std::endl;
  
  for (size_t i = 0; i < trees.size(); ++i) {
  	myfile << ale->p(trees[i]) << "\t" << ale2->p(trees[i]) << std::endl;
	cout  <<ale->p(trees[i]) << "\t" << ale2->p(trees[i]) << endl;
  }
  myfile.flush();
  myfile.close();
  delete ale;
  delete ale2;
}
  return 1;

}
