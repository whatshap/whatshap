#include "DP_matrix.h"

using namespace std;

/*
Computes the mutation costs according to the mut-function of the theoretical framework
*/
float mut(int refA, int refB, int hA, int hB, float c) {
	return c*(abs((refA-hA))+abs((refB-hB)));
}

/*
	performs the actual computation of the scoring matrix and returns this
*/
float** perform_DP(vector<int> H_A, vector<int> H_B, vector<vector<int>> Ref, float param_mut, unordered_set<int> E){
	int m = H_A.size();
	int n = Ref.size();

	//Initialize scoring matrix	
	float** Score = new float*[m];
	for(int i = 0; i < m; ++i)
		Score[i] = new float[n*n];
	//Initialize the arrays containing minima of columns
	float* minrightlist = new float[n];
	float* minleftlist = new float[n];
	for (int k = 0; k < n; ++k) {
		minrightlist[k] = INT_MAX;	
		minleftlist[k] = INT_MAX;	
	}
	
	int x = Ref.size();
	int y = Ref[0].size();
	
	//initialize the reference matrix
	int** ref = new int*[x];
	for(int i = 0; i < x; ++i)
		ref[i] = new int[y];
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
		 	ref[i][j] = Ref[i][j];	
		}
	}
	
	//converts the haplotype vectors to int arrays
	int* ha = new int[H_A.size()];
	int* hb = new int[H_B.size()];
	for (unsigned int i = 0; i < H_A.size(); i++) {
		ha[i] = H_A[i];	
		hb[i] = H_B[i];
	}
		
	float mingeneral = INT_MAX;
	//initialize first row of scoring matrix
	for(int left = 0; left < n; left++) {
		for (int right = 0; right < n; right++)  {
			int j = left*n + right;
			Score[0][j] = mut(Ref[left][0],Ref[right][0],H_A[0],H_B[0],param_mut);
			if (Score[0][left*n + right] < minrightlist[right]){
				minrightlist[right] = Score[0][left*n + right];
			}
			if (Score[0][left*n + right] < minleftlist[left]){
				minleftlist[left] = Score[0][left*n + right];
			}		
			if (Score[0][j] < mingeneral) mingeneral = Score[0][left*n + right];
		}
	}
	
	
	float  mutcost, switchcost = 0;
	float c00, c01, c10,c11,c02,c12,c21 = 0;
	int j = 0;
	for (int i = 1; i < m; ++i) {

		
		//fills the scoring matrix using the efficient computation of  O(m^2 x n)
		for (int left = 0; left < n; left++) {
			for (int right = 0; right < n; right++)  {
				j = left*n + right;
				mutcost = mut(ref[left][i],ref[right][i],ha[i],hb[i], param_mut);
				c00 = Score[i-1][j];
				c01 = 1+minrightlist[right];
				c10 = 1+minleftlist[left];
				c11 = 2+mingeneral;
				if (E.find(i-1) != E.end()) {
					c02 = Score[i-1][right*n + left];
					c12 = 1+minrightlist[left];
					c21 = 1+minleftlist[right];
					float mins[] = {c00,c01,c10,c11,c02,c12,c21};			
					switchcost = *min_element(mins,mins+7);
				}
				else {
					float mins[] = {c00,c01,c10,c11};
					switchcost = *min_element(mins, mins+4);				
				}
				Score[i][j] = switchcost + mutcost;					
			}
		}
		for (int k = 0; k < n; ++k) {
					minrightlist[k] = INT_MAX;	
					minleftlist[k] = INT_MAX;	
		}
		mingeneral = INT_MAX;
		
		//updating of the lists containing the column minima for later steps
		int j =0;
		for (int left = 0; left < n; left++) {
			for (int right = 0; right < n; right++)  {
				j = left*n + right;
				if (Score[i][j] < minrightlist[right]) minrightlist[right] = Score[i][j];				
				if (Score[i][j] < minleftlist[left]) minleftlist[left] = Score[i][j];				
				if (Score[i][j] < mingeneral) mingeneral = Score[i][j];
			}
		}
	}
	delete [] minrightlist;
	delete [] minleftlist;
	
	return(Score);
	for(int i = 0; i < m; ++i) {
    delete [] Score[i];
	}
	delete [] Score;
}	

tuple<float, int> minimum(float* current_list, int l, int r, bool is_terminal, int N) {
	/*
	Needed for the backtracing. Performs for a given column in current_list the backward computation
	 of the scoring function to determine the optimal value it origins from. 
	 Output: The minimum value in the preceding column and the according index
	*/
	vector<float> new_list;
	int size = N*N;
	int n = (int)sqrt(size);
	int switchval = 0;
	//backtracing from one column in the DP matrix to a column that was at the end of a block
	if (is_terminal) {
		for (int i = 0; i < size; i++) {			
			int l_ = (int)(i/n);
			int r_ = (int)(i%n);
			if ((l_ == r && r_==l) || (l_ == l && r_ == r))
				switchval = 0;
			else if ((l_ == r && r_ != l) || (l_ != r && r_ == l))
				switchval = 1;
			else if ((l_ == l && r_ != r) || (l_ != l && r_ == r))
				switchval = 1;
			else
				switchval = 2;
			new_list.push_back(current_list[i]+switchval);
		}
		//find optimal value when considering switch costs between columns
		float newmin = *min_element(new_list.begin(), new_list.end());
		int newmin_index = find(new_list.begin(), new_list.end(), newmin) - new_list.begin();
		return make_tuple(newmin, newmin_index);
	}
	//backtracing to a column that was not the end of a block
	else {
		for (int i = 0; i < size; i++) {
			int l_ = (int)(i/n);
			int r_ = (int)(i%n);
			if (l_ != l && r_ !=r)
				switchval = 2;
			else if (r_ != r && l_ == l)
				switchval = 1;
			else if (r_ == r && l_ != l)
				switchval = 1;
			else if (r_ == r && l_ == l) switchval = 0;
			new_list.push_back(current_list[i]+switchval);
		}
		int newmin = *min_element(new_list.begin(), new_list.end());
		int newmin_index = find(new_list.begin(), new_list.end(), newmin) - new_list.begin();
		return make_tuple(newmin, newmin_index);
	}		
}

/*
	Reads the input and calls the function for computation of the scoring matrix
*/
void compute_scoring(char* haplofile, char* E_file, char* panelfile, float mutation, char* out_costfile, char* out_pathfile){	
	//haplofile: A file with the haplotype strings extracted from the target vcf
	
	ifstream file(haplofile);
 	string haplo1;
 	string haplo2;
	getline(file, haplo1);
   getline(file, haplo2);
   
   //E_file: A file containing the indices for the block ending boundaries
  /* int counter =-1;
 		
 	ifstream infile(E_file);
  	string line;
	while(getline(infile,line)){
   	stringstream linestream(line);
    	string value;
   	while(getline(linestream,value,',')){
   		counter += 1;
		}
	}
	*/
	vector<int> vectorE;
	ifstream infile2(E_file);
	int counter2 = -1;
  	string line2;
	while(getline(infile2,line2)){
   	stringstream linestream(line2);
    	string value;
   	while(getline(linestream,value,',')){
   		counter2 += 1;
      	vectorE.push_back(stoi(value));
		}
	}
	unordered_set<int> e{begin(vectorE), end(vectorE)};
	
	//push strings for the haplotypes into vectors of ints
	vector<int> H_A;
	for (unsigned int i = 0; i < haplo1.size(); i++) {
		H_A.push_back(haplo1[i] - '0');	
	}
	vector<int> H_B;
	for (unsigned int i = 0; i < haplo2.size(); i++) {
		H_B.push_back(haplo2[i] - '0');	
	}
	
	//panelfile: The txt file containing the formerly created reference panel
	ifstream ref_panel(panelfile);
	string newline;
	//create a reference matrix
	vector<vector<int> > Ref;
	while (getline(ref_panel, newline)) {
		istringstream iss(newline);
		char ch;
		vector<int> row;
		while(iss >> ch) {
			row.push_back(ch-'0');		
		}
		Ref.push_back(row);
	}

	//read in the mutation parameter
		float param_mut = mutation;
	
	//perform the actual computation of the scoring matrix and return the matrix 'Score'	
	float** Score = perform_DP(H_A,H_B,Ref, param_mut,e);

	
	int end = H_A.size() -1;
	int n = Ref.size();

	//determines the starting point for the backtracing
	vector<int> path1;
	vector<int> path2;
	float mini = Score[end][0];
	int minarg = 0;
	int l = 0;
	int r = 0;
	vector<float> Score_vec;
	for (int i = 0; i < n*n; i++) {
		Score_vec.push_back(Score[end][i]);
		if (Score[end][i] < mini) {
			minarg = i;		
			mini = Score[end][i];
		}
	}

	//costfile: An output file where costs are written into	
	ofstream costfile;
	costfile.open(out_costfile);
	costfile << mini;
	
/* //for testing purposes only
	int doublecounter = 0;
	for (int i = 0; i < n*n; i++) {
		if (Score[end][i] == mini)
			doublecounter += 1;
	}
*/	

	//backtracing to create the paths through the scoring matrix
	l = (int)minarg/n;
	r = (int)minarg%n;
	path1.push_back(l);
	path2.push_back(r);
	while (end > 0) {
		Score_vec.clear();
		for (int i = 0; i < n*n; i++) {
			Score_vec.push_back(Score[end-1][i]);
		}
		minarg = get<1>(minimum(Score[end-1], l, r, e.find(end-1)!= e.end(), n) );
		l = (int)minarg/n;
		r = (int)minarg%n;

		path1.push_back(l);
		path2.push_back(r);
		end -= 1;
	}
	
	//paths were assembled from back, so they need to be reversed
	reverse(path1.begin(),path1.end());
	reverse(path2.begin(),path2.end());

	//outfile for the paths 	
	ofstream myfile;
   myfile.open(out_pathfile);
	for (auto &i : path1) {
		myfile << i << " ";
	}
	myfile << '\n';
	for (auto &i : path2) {
		myfile << i << ' ';
	}
	myfile << '\n';
}