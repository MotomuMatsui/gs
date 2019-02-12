/******************************************\
| Graph Splitting Method v2.4 (2019/02/12) |
|                                          |
| Copyright (c) 2015-2019 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2019)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

using namespace std;

void original(ifstream& ifs, unordered_map<string, int>& ORI, int const& size){

  string newick = "";

  string line;
  while(getline(ifs, line)){
    newick += line;
  }

  int N = newick.size();
  int B = 0; //brachet
  int D = 0;

  string id = "";

  unordered_map<int, vector<int>> clade;
  stringstream ss;

  for(int n=0; n<N; n++){
    auto p = newick[n];

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
      
      ORI[C] = 1;

      clade[B] = {};
      B--;
      D = 0;
    }
    else if(p == ','){
      if(id.size() > 0){
	for(int depth=1; depth<=B; depth++){
	  clade[depth].push_back(stoi(id));
	}
	id = "";
      }
      D = 0;
    }
    else if(p == ':'){
      D = 1;
    }
    else{
      if(D == 0){
	id += {p};
      }
    }
  }
}

void branch(string const& newick, unordered_map<string, int>& ORI, unordered_map<string, double>& ep, string& BRA, int const& ep_num, int const& size){

  BRA = "";
  int N = newick.size();
  int B = 0; //brachet
  int f = 1;
  string id = "";

  unordered_map<int, vector<int>> clade;
  stringstream ss;
  stringstream ss_BR;

  for(int n=0; n<N; n++){
    auto p = newick[n];

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
	int judge = (ORI[C]==1)? 1: 0;
	ss_BR << ep[C]/ep_num << "\t" << judge << "\n";	
      }
      else if(B==2 && f==1){
        int count  = 0;
        size_t pos = 0;
        while((pos = C.find("|", pos)) != string::npos){
          count ++;
          pos++;
        }

        if(count < size-1){
	  int judge = (ORI[C]==1)? 1: 0;
	  ss_BR << ep[C]/ep_num << "\t" << judge << "\n";	
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

  BRA = ss_BR.str();
}
