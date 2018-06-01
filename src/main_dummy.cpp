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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>
#include <unistd.h>
#include <random>   
#include <functional>

/*
/// format.cpp (File/Data-handling)
extern void readFASTA(ifstream&, ofstream&, ofstream&, int&);
extern void bl2mat(ifstream&, double*&, int const&);
extern void sc2nwk(int* const&, string&, int const&);
extern void addEP(string const&, string&, unordered_map<string, double>&, int const&, int const&);

/// mmseqs.cpp (Wrapper function of MMseqs)
extern void mmseqs(string const&, string const&, string const&, string const&);

/// gs.cpp (Core functions of GS method)
extern void GS(double* const&, int*&, int const&);
extern void EP(double* const&, unordered_map<string, double>&, function<double()>&, int const&);
*/

/// messages.cpp
extern void print_banner();
extern void print_usage(char*&);

//eigen.cpp (Wrapper function of LAPACKE/CBLAS package)
extern int eigen_lapack(double*&, double*&, int);

using namespace std;

int main(int argc, char* argv[]){
  
  /*/ Check the mmseqs command /*/
  
  //auto status = system("which mmseqs &> /dev/null");
  //if(WEXITSTATUS(status) != 0){
  //  /*PRINT*/ print_banner();
  //  /*PRINT*/ cerr << "mmseqs not found!\nCheck our web page (https://github.com/MotomuMatsui/gs) for more information" << endl;
  //  return -1;
  //}
  
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
    }
    else if(opt == 'r'){ // OK! (./gs -r 12345 IN.fst)
      if(regex_match(optarg, renum)){
	seed = atoi(optarg);
      }
    }
    else if(opt == 't'){ // OK! (./gs -t 4 IN.fst)
      if(regex_match(optarg, renum)){
	auto th = atoi(optarg);
	if(th >0){
	  threads = string(optarg);
	}
      }
    }
    else if(opt == 'm'){ // OK! (./gs -m 7.5 IN.fst)
      if(regex_match(optarg, renum)){
	auto sen = atof(optarg);
	if(1<=sen && sen<=7.5){
	  sensitivity = string(optarg);
	}
      }
    }
    else if(opt == 'h'){ // HELP message (./gs -h)
      return 0;
    }
    else if(opt == 'v'){ // Version (./gs -v)
      return 0;
    }
    else if(opt == 's'){ // SILENT mode (./gs -s -e 100 IN.fst)
      silence = 1;
    }
    else if (opt == '?'){
      if(optopt == 'e'){ // NG! (./gs IN.fst -e)
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
    return -1;
  }

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

  /*PRINT*/ print_banner();
  /*PRINT*/ print_usage(argv[0]);


  // Cluster size
  int N = 20;
  double* A = new double[N*N]();  
  double* z = new double[N*N]();     //z: eigenvector

  auto info = eigen_lapack(A, z, N); //LAPACK!

  cout << "Hello World!!" << endl;

  return 0;
}
