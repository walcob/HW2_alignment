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
  std::cout<<'\n';
  
  std::ofstream matrix;
  matrix.open("mattest.txt");
  
  int Nacross, Ndown; //number of collumns and Ndown in matrix
  std::string seq1, seq2;
  
  //read the matrix
  std::ifstream read, read_s1, read_s2;
  read_s1.open(argv[1]);
  read_s2.open(argv[2]);
  if(!read_s1.is_open()){
    std::cerr<<"File "<<argv[1]<<" could not be opened!\n";
    return -1;
  }
  else{
    seq1 = read_sequence(read_s1);
  }
  
  if(!read_s2.is_open()){
    std::cerr<<"File "<<argv[2]<<" could not be opened!\n";
    return -1;
  }
  else{
    seq2 = read_sequence(read_s2);
  }
   
  read.open(argv[3]);
  if(!read.is_open()){
    std::cerr<<"File "<<argv[3]<<" could not be opened!\n";
    return -1;
  }
  
  read>>Nacross>>Ndown;
  if(Nacross != seq1.size() || Ndown != seq2.size()){
      std::cerr<<"Alignment score files not the right dimensions.\n";
      return -1;
  }
  int S[Ndown][Nacross]; //score matrix
  //read scores
  matrix<<Ndown<<" "<<Nacross<<std::endl;
  std::string tmp;
  std::getline(read, tmp);
  matrix<<"Score Matrix: \n";
  for(int i = 0; i<Ndown; i++){
    for(int j = 0; j<Nacross; j++){
      std::getline(read, tmp);
      S[i][j]  = std::atoi(tmp.c_str());
      if(S[i][j]>=0) matrix<<" ";
      matrix<<S[i][j]<<" ";
    }
      matrix<<std::endl;
  }
  

  int A[Ndown][Nacross];  //alignment matrix

    
  /*for(int i = 1; i<=Nacross; i++)
      A[i][0]=seq1[i-1];
  for(int i = 1; i<=Ndown; i++)
      A[0][i]=seq2[i-1];    */
 
 
  int T[Ndown][Nacross];  //traceback matrix
  
  
  
  //Filling in gap
  A[0][0] = T[0][0] = 0;
  for(int i = 1; i < Nacross; i++){
      A[0][i]=0;
      T[0][i]=3;
  }
  for(int i = 1; i < Ndown; i++){
      A[i][0] = 0;
      T[i][0] = 2;
  }
  
  //meat of the program starts here.  
  //Recalculate adjusted alignment 
  //scores and fill in T
  for(int i = 1; i<Ndown-1; i++){
      for(int j = 1; j<Nacross-1; j++){
          A[i][j] = max(A[i-1][j-1]+S[i][j],A[i-1][j]-gap, A[i][j-1]-gap); 
          T[i][j] = argmax(A[i-1][j-1]+S[i][j],A[i-1][j]-gap, A[i][j-1]-gap);
      }
  }
	/*//last row
	for(int i = 1; i<Ndown; i++){
		A[i][Nacross] = max(A[i-1][Nacross-1]+S[i][Nacross], A[i-1][Nacross], A[i][Nacross-1]-gap);
		T[i][Nacross] = argmax(A[i-1][Nacross-1]+S[i][Nacross], A[i-1][Nacross], A[i][Nacross-1]-gap);
  }
  
  //last collumn
  for(int j = 1; j<Nacross; j++){
    A[Ndown][j] = max(A[Ndown-1][j-1]+S[Ndown][j], A[Ndown][j-1]-gap, A[Ndown-1][j]);
    T[Ndown][j] = argmax(A[Ndown-1][j-1]+S[Ndown][j], A[Ndown][j-1]-gap, A[Ndown-1][j]);
  }
  
  //last box
  A[Ndown][Nacross] = max(A[Ndown-1][Nacross-1] + S[Ndown][Nacross], A[Ndown][Nacross-1], A[Ndown-1][Nacross]);
  T[Ndown][Nacross] = argmax(A[Ndown-1][Nacross-1] + S[Ndown][Nacross], A[Ndown][Nacross-1], A[Ndown-1][Nacross]);
  */
  //test A
  matrix<<'\n'<<"Alignment Matrix: \n";
  for(int y = 0; y <Ndown; y++){
    for(int z = 0; z <Nacross; z++){
      if(A[y][z] >= 0) matrix<<" ";
      matrix<<A[y][z]<<" ";
    }
    matrix<<std::endl;
  }
  
  
  //test T
  
  matrix<<'\n'<<"Traceback Matrix: \n";
  for(int y = 0; y <Ndown; y++){
    for(int z = 0; z <Nacross; z++){
      matrix<<T[y][z]<<" ";
    }
    matrix<<std::endl;
  }
 /*
  //Traceback
  std::list <int> traceback;
  int i, j;
  i = Nacross-1;
  j = Ndown-1;
  while(i != 0 && j != 0){
    int k = T[i][j];
    traceback.push_back(k);
    if(k == 1){
      i--;
      j--;
    }
    else if (k == 2){
      i--;
    }
    else{
      j--;
    }
  }
  
  std::list <int>::iterator pos = traceback.end();
  pos--;
  std::string align1, align2;
  i = 0;
  while(i < seq1.size()){
    if(*pos == 1 || *pos == 2){
      align1 += seq1[i];
      i++;
    }
    else{
      align1 += '~';
    }
    pos--;
  }
  
  pos = traceback.end();
  pos--;
  
  j = 0;
  while(j < seq2.size()){
    if(*pos == 1 || *pos == 3){
      align2 += seq2[j];
      j++;
    }
    else{
      align2 += '~';
    }
    pos--;
  }
  
  std::ofstream result;
  result.open(argv[4]);
  if(!result.is_open()){
    std::cerr<<"File "<<argv[4]<<" could not be opened.\n";
    return -1;
  }
  result<<align1<<std::endl<<align2<<std::endl;*/
  return 0;
}