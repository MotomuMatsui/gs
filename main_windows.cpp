/******************************************\
| Graph Splitting Method v2.0 (2018/06/01) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <algorithm>
#include <cmath> 
#include <functional>
#include <fstream>
#include <iostream>
#include <iterator>
#include <lapacke.h>
#include <map>
#include <random>   
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <vector>

using namespace std;

/// format.cpp (File/Data-handling)
void readFASTA(ifstream&, ofstream&, ofstream&, int&);
void bl2mat(ifstream&, double*&, int const&);
void sc2nwk(int* const&, string&, int const&);
void addEP(string const&, string&, unordered_map<string, double>&, int const&, int const&);

/// mmseqs.cpp (Wrapper function of MMseqs)
void mmseqs(string const&, string const&, string const&, string const&);

/// gs.cpp (Core functions of GS method)
void GS(double* const&, int*&, int const&);
void EP(double* const&, unordered_map<string, double>&, function<double()>&, int const&);

/// messages.cpp
void print_banner();
void print_usage(char*&);

// sc.cpp (spectral clustering)
tuple<vector<int>,vector<int>> spectral_clustering(double* const&, vector<int> const&, int const&);

// gs_functions.cpp
void sedMATRIX(int*&, vector<int>&, vector<int> const&, int const, int const);
void sedVECTOR(vector<int>&, vector<int> const&, int const);
double simI(double* const&, vector<int> const&, int const);
int whichMIN(vector<double> const&);
double gev(double const&, double const&);

// sc_functions.cpp (Split one cluster into two subclusters)
void subMATRIX(double* const&, double*&, double*&, vector<int> const&, int const);
void whichCUT(double* const&, int const, double* const&, int*&, int&, int const);
void splitVECTOR(vector<int> const&, int const, int* const&, int const, vector<int>&, vector<int>&, int const);

// eigen.cpp (Wrapper function of LAPACKE/CBLAS package)
int eigen_lapack(double*&, double*&, int);

int main(int argc, char* argv[]){
  
  /*/ Check the mmseqs command /*/
  auto status = system("which mmseqs &> /dev/null");
  if(WEXITSTATUS(status) != 0){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "mmseqs not found!\nCheck our web page (https://github.com/MotomuMatsui/gs) for more information" << endl;
    return -1;
  }

  /*/ Getopt /*/
  int silence = 0;
  int ep_num  = 0;
  int seed    = 0;
  string threads = "1";
  string sensitivity = "7.5";
  int opt;
  regex renum(R"(^[\d\.]+$)"); // -e/-r/-t option requires an integer number
  opterr = 0;              // default error messages -> OFF

  while ((opt = getopt(argc, argv, "shve:r:t:m:")) != -1){
    if(opt == 'e'){ // OK! (./gs -e 100 IN.fst)
      if(regex_match(optarg, renum)){
	ep_num = atoi(optarg);
      }
      else{ // NG! (./gs -e hundred IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'r'){ // OK! (./gs -r 12345 IN.fst)
      if(regex_match(optarg, renum)){
	seed = atoi(optarg);
      }
      else{ // NG! (./gs -r one_two_three IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 't'){ // OK! (./gs -t 4 IN.fst)
      if(regex_match(optarg, renum)){
	auto th = atoi(optarg);
	if(th >0){
	  threads = string(optarg);
	}
	else{ // NG! (./gs -t four IN.fst)
	  /*PRINT*/ print_banner();
	  /*PRINT*/ cerr << "Option -t requires an integer argument (>=1).\n" << endl;
	  /*PRINT*/ print_usage(argv[0]);
	  return -1;
	}
      }
      else{ // NG! (./gs -t four IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -t requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'm'){ // OK! (./gs -m 7.5 IN.fst)
      if(regex_match(optarg, renum)){
	auto sen = atof(optarg);
	if(1<=sen && sen<=7.5){
	  sensitivity = string(optarg);
	}
	else{ // NG! (./gs -m 10 IN.fst)
	  /*PRINT*/ print_banner();
	  /*PRINT*/ cerr << "Option -m requires a double number argument [1, 7.5].\n" << endl;
	  /*PRINT*/ print_usage(argv[0]);
	  return -1;
	}
      }
      else{ // NG! (./gs -m seven IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -m requires a double number argument [1, 7.5].\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'h'){ // HELP message (./gs -h)
      /*PRINT*/ print_banner();
      /*PRINT*/ print_usage(argv[0]);
      return 0;
    }
    else if(opt == 'v'){ // Version (./gs -v)
      /*PRINT*/ print_banner();
      return 0;
    }
    else if(opt == 's'){ // SILENT mode (./gs -s -e 100 IN.fst)
      silence = 1;
    }
    else if (opt == '?'){
      if(optopt == 'e'){ // NG! (./gs IN.fst -e)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
      else{ // NG! (./gs -Z)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << argv[0] << ": invalid option\n" <<  endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
  }

  /*/ Input file name /*/
  string input = "";
  if(optind < argc){ // OK! (./gs -e 100 IN.fst)
    input = argv[optind];
  }
  else{ // NG! (./gs -e 100)
    /*PRINT*/ cerr << argv[0] << " requires an input file (fasta format).\n" << endl;
    /*PRINT*/ print_banner();
    /*PRINT*/ print_usage(argv[0]);
    return -1;
  }
  
  /*/ Variables /*/
  double* W;     // Sequence similarity matrix
  int size;      // Row size of W (W is a synmetry matrix)
  int* gs;       // Result of GS method
  string newick; // GS tree (without EP values)

  /*/ File I/O /*/
  regex re(R"(\.[^\.]+$)"); // Extention of the input file
  auto original_fasta = string(input);
  auto annotation_txt = regex_replace(original_fasta, re, "_annotation.txt");
  auto simple_fasta   = regex_replace(original_fasta, re, "_simple.fst");
  auto mmseqs_result  = regex_replace(original_fasta, re, "_mmseqs.txt");
  
  ifstream ifs1(original_fasta); // Fasta file (original)
  ofstream ofs1(annotation_txt); // Annotation file
  ofstream ofs2(simple_fasta);   // Fasta file (simple)
  
  if(ifs1.fail()){
    /*PRINT*/ cerr << "\nCannot access " << original_fasta << "!" << endl;
    return -1;
  }
  if(ofs1.fail()){
    /*PRINT*/ cerr << "\nCannot create " << annotation_txt << "!" << endl;
    return -1;
  }
  if(ofs2.fail()){
    /*PRINT*/ cerr << "\nCannot create " << simple_fasta << "!" << endl;    
    return -1;
  }

  /*/ Parsing fasta file /*/
  readFASTA(ifs1, ofs1, ofs2, size); 
    // ifs1: INPUT (original fasta file)
    // ofs1: OUTPUT (annotation file)
    // ofs2: OUTPUT (simplified fasta file)
    // size: # of sequence = row size of sequence similarity matrix
  
  /*/ Parameters /*/  
  if(!silence){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Settings:" << endl;
    /*PRINT*/ cerr << "-Input" << endl;
    /*PRINT*/ cerr << "  file = " << input << endl;
    /*PRINT*/ cerr << "  # of sequences = " << size << endl << endl;

    /*PRINT*/ cerr << "-MMseqs" << endl;
    /*PRINT*/ cerr << "  Sentitivity = " << sensitivity << endl;
    /*PRINT*/ cerr << "  # of threads = " << threads << endl << endl;

    /*PRINT*/ cerr << "-EP method" << endl;
    if(seed>0){
      /*PRINT*/ cerr << "  Random seed = " << seed << endl;
    }
    else{
      /*PRINT*/ cerr << "  Random seed = " << "a random number (default)" << endl;
    }
    /*PRINT*/ cerr << "  # of iterations = " << ep_num << endl << endl;

    /*PRINT*/ cerr << "Progress:" << endl;
  }
  
  /*/ Executing MMSeqs /*/
  /*PRINT*/ if(!silence) cerr << "-MMseqs\n  " << size << "x" << size << " pairwise alignment\n" << "  searching...\r" << flush;
  mmseqs(simple_fasta, mmseqs_result, threads, sensitivity);
    // simple_fasta: INPUT (multiple fasta file) 
    // mmseqs_result: OUTOUT (result file of MMseqs)
    // threads: parameter (threads number for MMseqs)
    // sensitivity: parameter (sensitivity for MMseqs)

  ifstream ifs2(mmseqs_result); // MMseqs result file
  if(ifs2.fail()){
    /*PRINT*/ cerr << "  Cannot access " << mmseqs_result << "!" << endl;
    return -1;
  }
  else{
    /*PRINT*/ if(!silence) cerr << "  done.        " << endl << endl;
  }

  /*/ Reading Data /*/
  bl2mat(ifs2, W, size);
    // ifs2: INPUT (result file of MMseqs)
    // W: OUTPUT (sequence similarity matrix)

  /*/ GS method (stepwise spectral clustering) /*/
  /*PRINT*/ if(!silence) cerr << "-GS method\n" << "  executing...\r" << flush;
  GS(W, gs, size);
    // W: INPUT (sequence similarity matrix)
    // gs: OUTPUT (result of stepwise spectral clustering)

  /*PRINT*/ if(!silence) cerr << "  done.         " << endl << endl;

  /*/ Generating GS tree Newick based on the spectral clustering /*/
  sc2nwk(gs, newick, size);
    // gs: INPUT (result of stepwise spectral clustering)
    // newick: OUTPUT (GS tree [newick format])

  /*/ EP method /*/
  if(ep_num>0){
    unordered_map<string, double> ep;
    string newick_EP; // GS+EP tree

    // Random number generator (Uniform distribution->Mersenne Twister)
    function<double()> R;
    uniform_real_distribution<double> urd(0,1);    // uniform distributed random number

    if(seed>0){
      mt19937 mt(static_cast<unsigned int>(seed)); // mersenne twister
      R = bind(urd, ref(mt));                      // random number generator    
    }
    else{
      random_device rd;                            // random seed
      mt19937 mt(rd());                            // mersenne twister
      R = bind(urd, ref(mt));                      // random number generator        
    }    

    /*PRINT*/ if(!silence) cerr << "-EP method" << endl;

    for(int n=1; n<=ep_num; n++){
      /*PRINT*/ if(!silence) cerr << "  " << n << "/" << ep_num << " iterations" << "\r"<< flush;
      EP(W, ep, R, size);
        // W: INPUT (sequence similarity matrix)
        // ep: OUTPUT (result of Edge Perturbation method)
        // R: random number generator
    }
    
    /*PRINT*/ if(!silence) cerr << "\n  done." << endl;
    /*PRINT*/ if(!silence) cerr << "\n------------------------------------------\n" << endl;

    addEP(newick, newick_EP, ep, ep_num, size);
      // newick: INPUT (GS tree [newick format])
      // newick_EP: OUTPUT (GS+EP tree [newick format])
      // ep: INPUT (result of Edge Perturbation method)
      // ep_num: INPUT (# of Edge Perturbation method)

    /*/ GS tree WITH EP values ->STDOUT /*/
    cout << newick_EP << endl;
  }
  else{ // skip the EP method
    /*PRINT*/ if(!silence) cerr << "\n------------------------------------------\n" << endl;

    /*/ GS tree WITHOUT EP values ->STDOUT /*/
    cout << newick << endl;
  }

  delete[] W;
  delete[] gs;

  return 0;
}

/*
 *
 * gs.cpp
 *
 */

void GS(double* const (&W), int* (&step), int const& size){

  vector<double> simL(size,3);
  vector<int> res(size,1);

  step = new int[size*size]();
  sedMATRIX(step, res, res, 1, size);

  int gK   = 1;
  int gMin = 1;
  vector<int> a;
  vector<int> b;

  while(gK < size){
    //Spectral clustering
    tie(a,b) = spectral_clustering(W, res, gMin);

    //Record the result of spectral clustering
    gK++;
    sedMATRIX(step, res, b, gK, size);
    
    //Decision of the cluster which will be analyzed at the next step
    simL[gMin-1] = simI(W, res, gMin); //cluster a
    simL[gK-1]   = simI(W, res, gK);   //cluster b

    gMin = whichMIN(simL);
  }
}

/*
 *
 * ep.cpp
 *
 */

void EP(double* const (&W), unordered_map<string, double>& ep, function<double()>& R, int const& size){

  // Edge perturbation
  int N = size*size;
  double* E = new double[N]();
  for(int n = 0; n < N; n++){
    double a = W[n];
    if(a==0){
      E[n] = 0; // Graph topology is not changed
    }
    else{
      auto b = gev(R(), a);
      E[n] = 
	(b>1)? 1: // b: Similarity scores perturved according to Generalized Extreme Value distribution
	(b<0)? 0: // E[n]: Perturved sequence similarity score
	       b; // 0 <= E[n] <= 1
    }
  }

  //GS method + EP method
  vector<double> simL(size,3);
  vector<int> res(size,1);

  int gK   = 1;
  int gMin = 1;
  vector<int> a;
  vector<int> b;

  stringstream ss;

  while(gK < size){
    //Spectral clustering
    tie(a,b) = spectral_clustering(E, res, gMin);
    
    ss.str("");
    sort(a.begin(), a.end());
    for(int n : a){
      ss << n+1 << "|";
    }
    ep[ss.str()]++;

    ss.str("");
    sort(b.begin(), b.end());
    for(int n : b){
      ss << n+1 << "|";
    }
    ep[ss.str()]++;

    //Record the result of spectral clustering
    gK++;
    sedVECTOR(res, b, gK);
    
    //Decision of the cluster which will be analyzed at the next step
    simL[gMin-1] = simI(E, res, gMin); //cluster a
    simL[gK-1]   = simI(E, res, gK);   //cluster b

    gMin = whichMIN(simL);
  }

  //Free
  delete[] E;
}

/*
 *
 * sc.cpp
 *
 */

tuple<vector<int>,vector<int>> spectral_clustering(double* const& oW, vector<int> const& res, int const& num){

  // Cluster size
  int N = (int)count(res.begin(), res.end(), num);
  
  double* W = new double[N*N]();
  double* D = new double[N]();
  double* A = new double[N*N]();  

  //Rayleigh quotient of the submatrix
  subMATRIX(oW, W, D, res, num);

  //A = E-D*W*D
  for(int x = 0; x < N; x++){
    for(int y = 0; y < x; y++){
      auto p = x*N+y;
      A[p] = -D[x]*D[y]*W[p];
    }
    auto p = x*N+x;
    A[p] = 1-D[x]*D[x]*W[p];
  }

  //Eigenvalue analysis
  double* z = new double[N*N]();     //z: eigenvector
  auto info = eigen_lapack(A, z, N); //LAPACK!
  if(info != 0){}                    //Error occurs if info>0

  //Split submatrix into two parts according to the result of spectral clustering
  int* qi = new int[N]();
  int cut = 0;
  vector<int> a, b;
  whichCUT(z, 0, D, qi, cut, N);
  splitVECTOR(res, num, qi, cut, a, b, N);
  
  //free
  delete[] A;
  delete[] W;
  delete[] D;
  delete[] z;
  delete[] qi;

  return forward_as_tuple(a,b);
}

/*
 *
 * gs_functions.cpp
 *
 */

void sedMATRIX(int* (&step), vector<int>& vec, vector<int> const& pos, int const num, int const N){
  for(auto p: pos){
    vec[p] = num;
  }
  for(int n = 0; n<N; n++){
    step[n*N+num-1] = vec[n];
  }
}

void sedVECTOR(vector<int>& vec, vector<int> const& pos, int const num){
  for(auto p: pos){
    vec[p] = num;
  }
}

//simI function
double simI(double* const (&W), vector<int> const& res, int const num){
  int N    = res.size();
  int s    = count(res.begin(), res.end(), num);
  double D = 0;

  if(s<=1){
    D = 2;
  }
  else{
    for(int ro=0; ro < N; ro++){
      if(res[ro] == num){
        for(int co=0; co < ro; co++){
          if(res[co] == num){
            D += 2*W[ro*N+co];
          }
        }
      }
    }

    D /= (s*(s-1));
  }

  return D;
}

//Minimum factor
int whichMIN(vector<double> const& vec){
  int N = vec.size();
  int p = 0;
  double old = 100;
  for(int n = 0; n < N; n++){
    p   = (vec[n]<old)? n: p;
    old = vec[n];
  }

  return(p+1);
}

//Generalized Extreme Value function (inverse function)                                                                        
double gev(double const& x, double const& mu){
  double theta = mu*(1-mu)/3;
  double gamma = exp(-3*mu)-1;

  if(gamma == 0){ //Gummbel distribution                                                                                        
    return mu - theta*log(-log(x));
  }
  else{
    return mu +( pow(-log(x),-gamma)-1 )*theta/gamma;
  }
}

/*
 *
 * sc_functions.cpp
 *
 */

void subMATRIX(double* const (&Wo), double* (&W), double* (&D), vector<int> const& res, int const num){

  int N = res.size();
  int s = count(res.begin(), res.end(), num);

  int rn   = 0;
  int cn   = 0;
  double p;
  for(int ro=0; ro < N; ro++){
    if(res[ro] == num){
      cn=0;
      for(int co=0; co < ro; co++){
        if(res[co] == num){
          p = Wo[ro*N+co];
          W[rn*s+cn] = p;
          W[cn*s+rn] = p;
          D[rn]     += p;
          D[cn]     += p;
          cn++;
        }
      }
      p = Wo[ro*N+ro];
      W[rn*s+rn] = p;
      D[rn]     += p;

      rn++;
    }
  }

  for(int p = 0; p<s; p++){
    D[p] = 1/sqrt(D[p]);
  }
}

//Best cut position
typedef pair<double, int> P;
bool comp(const P &a, const P &b){
  return a.first > b.first;
}

void whichCUT(double* const (&z), int const col, double* const (&D), int*& qi, int& cut, int const N){
  vector<P> pairs(N);
  double* qs = new double[N]();

  //Inner product: <eigenvector,D>
  for(int p = 0; p < N; p++){
    pairs[p] = make_pair(z[col+p*N]*D[p], p);
  }

  //Sorting values
  sort(pairs.begin(), pairs.end(), comp);
  for(int p = 0; p < N; p++){
    qs[p] = pairs[p].first;
    qi[p] = pairs[p].second;
  }

  cut = 0;
  double best = 0;
  for(int n = 0; n < N-1; n++){
    auto def = abs(qs[n] - qs[n+1]);
    if(def > best){
      cut  = n;
      best = def;
    }
  }

  delete[] qs;
}

void splitVECTOR(vector<int> const& res, int const num, int* const (&qi), int const cut, vector<int>& a, vector<int>& b, int const N){

  vector<int> index;
  int n = 0;
  for(int v: res){
    if(v == num){
      index.push_back(n);
    }
    n ++;
  }

  for(int p = 0; p <=cut; p++){
    a.push_back(index[qi[p]]);
  }
  for(int p = cut+1; p < N; p++){
    b.push_back(index[qi[p]]);
  }
}

/*
 *
 * mmseqs.cpp
 *
 */

void mmseqs(string const& fst, string const& output, string const& threads, string const& sensitivity){
  //mmseqs command
  const string mmseqs = "mmseqs";

  // File I/O
  regex re(R"(_simple\.fst$)"); // Matching to the extention of the input file
  auto DB     = regex_replace(fst, re, "_DB");
  auto result = regex_replace(fst, re, "_result");
  auto tmp    = regex_replace(fst, re, "_tmp");

  //Commands
  auto cmd1 = mmseqs+" createdb "    +fst+" "+DB+" > /dev/null";
  auto cmd2 = mmseqs+" search "      +DB+" "+DB+" "+result+" "+tmp+" --threads "+threads+" -e 10 -s "+sensitivity+" > /dev/null";
  auto cmd3 = mmseqs+" convertalis " +DB+" "+DB+" "+result+" "+output+" --format-mode 0 > /dev/null";

  mkdir(tmp.c_str(), S_IRWXU | S_IRGRP | S_IROTH); //mkdir tmp
  auto info1 = system(cmd1.c_str()); if(info1>0){} //MMseqs: createdb
  auto info2 = system(cmd2.c_str()); if(info2>0){} //MMseqs: search
  auto info3 = system(cmd3.c_str()); if(info3>0){} //MMseqs: convertalis

  //Delete intermediate files
  string extension[6] = {
    "", ".dbtype", ".lookup",
    ".index", "_h.index", "_h"};
  
  for(auto e : extension){
    auto file1 = DB     + e;
    auto file2 = result + e;
    remove(file1.c_str());
    remove(file2.c_str());
  }
}

/*
 *
 * eigen.cpp
 *
 */

int eigen_lapack(double* (&A), double* (&z), int N){

  //Variables for LAPACKE_dsyevr function
  double  vl     = 1;
  double  vu     = N;
  int     il     = 2; // Graph Splitting method (Spectral clustering)
  int     iu     = 2; // requires only the 2nd smallest eigenvalue
  double  abstol = 0;

  int* m      = new int[iu-il+1];
  double* w   = new double[N];
  int* isuppz = new int[(iu-il+1)*N];
  
  //Eigenvalue analysis using LAPACK/BLAS (real symmetric matrix)
  auto info = LAPACKE_dsyevr(
			     LAPACK_ROW_MAJOR, // Matrix layout (Row major layout)
			     'V',              // Jobs (V: eigenvalues and eigenvectors)
			     'I',              // Range (I: the IL-th through IU-th eigenvalues will be found)
			     'L',              // Stored triangle matrix (L: Lower triangle of A is stored)
			     N,                // Order of A
			     A,                // Real symmetric matrix
			     N,                // Leading dimension of A
			     vl,vu,            // (These arguments are NOT referenced, but required as dummies)
			     il,iu,            // Index of the smallest (il) and largest (iu) eigenvalues to be returned
			     abstol,           // Absolute error tolerance
			     m,                // Total number of eigenvalues
			     w,                // Eivenvalues
			     z,                // Eigenvectors
			     N,                // Leading dimension of z
			     isuppz            // Support value of z
			     );
  
  if(info != 0){
    //Error occurs if info>0
  }

  //Free
  delete[] m;
  delete[] w;
  delete[] isuppz;

  return info;
}

/*
 *
 * format.cpp
 *
 */

// Parsing multiple fasta file
void readFASTA(ifstream& ifs, ofstream& ofs1, ofstream& ofs2, int& row){

  //Input & Output
  int id = 1;

  string line;
  while(getline(ifs, line)){
    if(line[0] == '>'){
      auto name = line.substr(1);
      ofs1 << id << "\t" << name << "\n";
      ofs2 << '>' << id << "\n";

      id ++;
    }
    else{
      ofs2 << line << "\n";
    }
  }

  // # of sequence (= row size of sequence similarity matrix)
  row = id - 1;

  ofs1.close();
  ofs2.close();
}

// Parsing mmseqs result file
void bl2mat(ifstream& ifs, double* (&W), int const& size){

  //File I/O
  double* S = new double[size*size](); // Bit score matrix
          W = new double[size*size](); // Sequence similarity matrix

  //Reading lines
  int x = 1; // Sequence IDs (query)
  int y = 1; // Sequence IDs (target)
  double b = 0; // Bit scores calculated by MMseqs
  string line, chr;
  while(getline(ifs, line)){
    //Split lines
    istringstream stream(line);

    int pos = 0;
    while(getline(stream, chr, '\t')){
      if(pos == 0){
	x = stoi(chr) -1;
      } 
      else if(pos == 1){
	y = stoi(chr) -1;
      }
      else if(pos == 11){
	b = stod(chr);
      }
      pos ++;
    }

    auto p = x*size+y;
    auto s = S[p];
    S[p] = s>b? s:b;
  }  

  for(int x = 0; x < size; x ++){
    for(int y = x+1; y < size; y ++){
      auto xy = S[x*size+y];
      auto yx = S[y*size+x];      
      auto xx = S[x*size+x];
      auto yy = S[y*size+y];
      
      auto comp = xy>yx?  xy:yx;
      auto self = xx>yy?  xx:yy;
      auto sss  = self>0? comp/self:0;

      W[x*size+y] = sss;
      W[y*size+x] = sss;
    }
    W[x*size+x] = 1;
  }
}

// Parsing result file generated by GS method
void sc2nwk(int* const& W, string& newick, int const& size){

  //Input & Output
  int* tr     = new int[size+1](); 
  int* change = new int[size*3]();

  //Reading lines
  for(int r=0; r<size; r++){
    int cur = 1;
    int old = 1;
    for(int c=0; c<size; c++){
      cur = W[r*size+c];
      if(cur != old){
	change[c*3]   = old; // previous membership
	change[c*3+1] = cur; // current membership
	change[c*3+2] ++;    // # of cluster member
      }
      old = cur;
    }
    tr[W[(r+1)*size-1]] = r+1;
  }

  int* nwk = new int[size+1](); for(int p=0; p<size; p++){nwk[p] = 1;}
  int* bra = new int[2*(size+1)](); // # of brachets
  int* pos = new int[2*(size+1)](); // left and right boundary of each cluster
  int* com = new int[size+1]();     // position of commas

  bra[0]        = 1;
  bra[size*2+1] = 1;
  pos[2]        = 0;
  pos[3]        = size;

  for(int p=1; p<size; p++){
    auto old = change[p*3];
    auto cur = change[p*3+1];
    auto num = change[p*3+2];

    auto ps  = pos[old*2];
    auto ls  = pos[old*2+1];
    
    pos[old*2]   += num; // 1 1 1 1 1 -> (2 2),(1 1 1) ... left boundary of cluster 1 = 2
    pos[old*2+1] -= num; // 1 1 1 1 1 -> (2 2),(1 1 1) ... # of cluster 1 = 3
    pos[cur*2]   = ps;   // 1 1 1 1 1 -> (2 2),(1 1 1) ... left boundary of cluster 2 = 0
    pos[cur*2+1] = num;  // 1 1 1 1 1 -> (2 2),(1 1 1) ... # of cluster 2 = 2
    com[ps+num]  = 1;    // 1 1 1 1 1 -> (2 2),(1 1 1) ... position of comma = 2
    
    for(int n=0; n<num; n++){
      nwk[ps+n] = cur;   // 1 1 1 1 1 -> (2 2),(1 1 1) ... replacing 1 -> 2
    }
    
    bra[ps*2]         ++; // (2 2
    bra[(ps+num)*2+1] ++; //     )
    bra[(ps+num)*2]   ++; //      (1 1 1
    bra[(ps+ ls)*2+1] ++; //            )
  }

  newick = "";
  for(int b=0; b<=size; b++){
    bra[b*2]   --; // (
    bra[b*2+1] --; // )

    string right = (bra[b*2+1]>0)? string(bra[b*2+1], ')') : "";
    string left  = (bra[b*2]  >0)? string(bra[b*2],   '(') : "";
    string comma = (com[b]    >0)? string(com[b],     ',') : "";
    string node  = (tr[nwk[b]]>0)? to_string(tr[nwk[b]])   : "";

    newick += right+comma+left+node;
  }
  newick += ";";

  //Free
  delete[] tr;
  delete[] change;
  delete[] nwk;
  delete[] bra;
  delete[] pos;
  delete[] com;
}

void addEP(string const& newick, string& newick_EP, unordered_map<string, double>& ep, int const& ep_num, int const& size){

  newick_EP = "";
  int N = newick.size();
  int B = 0; //brachet
  int f = 1;
  string id = "";

  unordered_map<int, vector<int>> clade;
  stringstream ss;
  stringstream ss_EP;

  for(int n=0; n<N; n++){
    auto p = newick[n];
    ss_EP << p;

    if(p == '('){
      B++;
    }
    else if(p == ')'){

      if(id.size() > 0){
	for(int depth=1; depth<=B; depth++){
	  clade[depth].push_back(stoi(id));
	}
	id = "";
      }

      ss.str("");
      sort(clade[B].begin(), clade[B].end());
      for(int n : clade[B]){
	ss << n << "|";	
      }

      string const C = ss.str();

      if(B>2){
	ss_EP << ep[C]/ep_num;
      }
      else if(B==2 && f==1){
	int count  = 0;
	size_t pos = 0;
	while((pos = C.find("|", pos)) != string::npos){
	  count ++;
	  pos++;
	}

	if(count < size-1){
	  ss_EP << ep[C]/ep_num;
	}

	f=0;
      }
 
      clade[B] = {};
      B--;
    }
    else if(p == ','){
      if(id.size() > 0){
	for(int depth=1; depth<=B; depth++){
	  clade[depth].push_back(stoi(id));
	}
	id = "";
      }
    }
    else{
      id += {p};
    }
  }

  newick_EP = ss_EP.str();
}

/*
 *
 * messages.cpp
 *
 */

void print_banner(){

  string banner = 
    "------------------------------------------\n"
    " Graph Splitting Method v2.0 (2018/06/01) \n"
    "                                          \n"
    "   Copyright (c) 2018 Motomu Matsui       \n"
    "   Systematic Biology, xx:xx-xxx, 2018    \n"
    "                                          \n"
    "   http://gs.bs.s.u-tokyo.ac.jp/          \n"
    "------------------------------------------\n";

  cerr << banner << endl;
}

void print_usage(char*& program){
  cerr << "Usage: " << program << " [-e INTEGER(>=0)] [-r INTEGER(>0)] [-t INTEGAR(>0)] [-s] [-h] [-v] input > output" << endl;
  cerr << "-e " << "the number of replicates for EP method. Default: 0" << endl;
  cerr << "-r " << "the random seed number for EP method. Default: random number" << endl;
  cerr << "-t " << "the number of threads for MMseqs. Default: 1" << endl;
  cerr << "-m " << "sensitivity for MMseqs. Default: 7.5" << endl;
  cerr << "-s " << "silent mode: do not report progress. Default: Off" << endl;
  cerr << "-h " << "show help messages. Default: Off" << endl;
  cerr << "-v " << "show the version. Default: Off" << endl;
}
