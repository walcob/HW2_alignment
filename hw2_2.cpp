#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <cstdlib>

int max(int a, int b, int c){
    if (a>b && a>c) return a;
    else if (b>a && b>c) return b;
    else return c;
}

int argmax(int a, int b, int c){
    if(max(a,b,c) == a) return 1;
    else if(max(a,b,c)== b) return 2;
    else return 3;
}
    
//this was borrowed directly from HW1
//take the FASTA sequence and convert it into a string with no blank spaces
std::string read_sequence(std::istream &seq){
    std::string aa="RHKDESTNQCGPAVILMFYW";
	std::string sequence;
	std::vector <std::string> lines;
	std::string line;
	while(!seq.eof()){
		std::getline(seq, line);
		lines.push_back(line);
	}
	
	//append each line into the sequence ignoring lines that begin with ">"
	int i = 0;
	while(i<lines.size()){
		if(lines[i][0]=='>') i++;
		else{
			sequence += lines[i];
			i++;
		}
	}
	
	//eliminate blank space characters from sequence
	std::string::iterator pos = sequence.begin();
	while(pos != sequence.end()){
		if(std::isspace(*pos)) sequence.erase(pos);
		else{
			//convert lowercase characters to capital letters
			*pos = toupper(*pos);
			pos++;
		}
        
	}
    
    //added this part to ensure that entire sequence is amino acid
    std::string::iterator a = sequence.begin();
    while(a!= sequence.end()){
        if(aa.find(*a) == std::string::npos){
            sequence.erase(i);
        }
        else{
            a++;
        }
    }
    
	return sequence;
}
int main(int argc, char* argv[]){ // ./hw2.exe seq1 seq2 score_matrix result
  
  int gap;
  std::cout<<"Enter gap penalty: ";
  std::cin>>gap;
  std::cout<<std::endl;
  
  std::ofstream matrix;
  matrix.open("mattest.txt");
  
  int Nacross, Ndown; //number of collumns and Ndown in matrix
  std::string seq1, seq2;
  
    //read the matrix
  std::ifstream read, read_s1, read_s2;
  read_s1.open(argv[1]);
  read_s2.open(argv[2]);
  if(!read_s1.is_open()){
    std::cerr<<"File "<<argv[1]<<" could not be opened!"<<std::endl;
    return -1;
  }
  else{
    seq1 = read_sequence(read_s1);
  }
  
  if(!read_s2.is_open()){
    std::cerr<<"File "<<argv[2]<<" could not be opened!"<<std::endl;
    return -1;
  }
  else{
    seq2 = read_sequence(read_s2);
  }
   
  read.open(argv[3]);
  if(!read.is_open()){
    std::cerr<<"File "<<argv[3]<<" could not be opened!"<<std::endl;
    return -1;
  }
  
  read>>Nacross>>Ndown;
  if(Nacross != seq1.size() || Ndown != seq2.size()){
      std::cerr<<"Alignment score files not the right dimensions."<<std::endl;
      return -1;
  }
  
  int S[Nacross+1][Ndown+1];
  std::string tmp;
  std::getline(read, tmp);
  for(int j = 1; j<=Ndown; j++){
    for(int i = 1; i<=Nacross; i++){
      std::getline(read, tmp);
      S[i][j] = std::atoi(tmp.c_str());
    }
  }
  
  int A[Nacross+1][Ndown+1];
  int T[Nacross+1][Ndown+1];
  T[0][0]=0;
  A[0][0]=0;
  
  //step 4
  for(int j = 1; j<=Ndown; j++){
    A[0][j] = 0;
    T[0][j] = 3;
  }
  
  for(int i = 1; i<=Nacross; i++){
    A[i][0]=0;
    T[i][0]=2;
  }
 
  //step 5
  for(int i = 1; i<Nacross; i++){
    for(int j = 1; j<Ndown; j++){
      A[i][j] = max(A[i-1][j-1]+S[i][j], A[i-1][j]-gap, A[i][j-1]-gap);
      T[i][j] = argmax(A[i-1][j-1]+S[i][j], A[i-1][j]-gap, A[i][j-1]-gap);
    }
  }
  
  //step 6, last row
  for(int i = 1; i<=Nacross; i++){
    A[i][Ndown] = max(A[i-1][Ndown-1]+S[i][Ndown], A[i-1][Ndown], A[i][Ndown-1]-gap);
    T[i][Ndown] = argmax(A[i-1][Ndown-1]+S[i][Ndown], A[i-1][Ndown], A[i][Ndown-1]-gap);
  }
  
  //step 7, last collumn
  for(int j = 1; j<=Ndown; j++){
    A[Nacross][j] = max(A[Nacross-1][j-1] + S[Nacross][j], A[Nacross-1][j]-gap, A[Nacross][j-1]);
    T[Nacross][j] = argmax(A[Nacross-1][j-1] + S[Nacross][j], A[Nacross-1][j]-gap, A[Nacross][j-1]);
  }  
  
  //step 8, Last box
  A[Nacross][Ndown] = max(A[Nacross-1][Ndown-1] + S[Nacross][Ndown], A[Nacross-1][Ndown], A[Nacross][Ndown-1]);
  T[Nacross][Ndown] = argmax(A[Nacross-1][Ndown-1] + S[Nacross][Ndown], A[Nacross-1][Ndown], A[Nacross][Ndown-1]);
  
  //step 9, Traceback
  std::list <int> traceback;
  int a,b;
  a = Nacross;
  b = Ndown;
  while( a >= 0 && b >= 0){
    traceback.push_back(T[a][b]);
    if(T[a][b] == 3)
      b--;
    else if(T[a][b]==2)
      a--;
    else{
      a--;
      b--;
    }
  }
  traceback.pop_back();

 std::list <int>::iterator pos = traceback.end();
 pos--;
 std::string align1, align2;
 int i, j;
 i = j = 0;
 int k =0;
 int l = max(seq1.size(), seq2.size(), 0);
 while(k<l){
  k++;
  if(*pos == 1){
    align1 += seq1[i];
    align2 += seq2[j];
    i++;
    j++;
    pos--;
  }
  else if (*pos == 2){
    align1 += seq1[i];
    align2 += '~';
    i++;
    pos--;
  }
  else{
    align1 += '~';
    align2 += seq2[j];
    j++;
    pos--;
  }
 }
  
  std::ofstream results;
  results.open(argv[4]);
  if(!results.is_open()){
    std::cerr<<"File "<<argv[4]<<" could not be opened."<<std::endl;
    return -1;
  }
  results<<align1<<std::endl<<align2;
  results.close();
 

  return 0;
}
