

/***************
 *Written by David Gfeller

 *The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

 *For any question or commercial use, please ask david.gfeller@unil.ch

 *Copyright (2022) David Gfeller
 ****************/

#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atoi */
#include <string.h>
#include <math.h> 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>

using namespace std;

int map(char c);
void find_pos(int *p, char* h, int l);

// This loads the predictions of MixMHCpred and gives as output the predictions of PRIME.
// To be even faster, we should include this directly in the C++ code of MixMHCpred.

int main(int argc, char ** argv){

    int Lmin=8;
    int Lmax=14;
    
    int N=20;
    
    string alphabet="ACDEFGHIKLMNPQRSTVWY";
    char letter[N];
    strcpy(letter, alphabet.c_str());
    
    char *lib_dir;
    int nh;
    int Np;
    int rd;
    int ct;
    string s;
    
    int inc;
    char **alleles;
    char **alleles_map;

    char *output_file;
    char *input_file;
    char *affinity_dir;
    
    lib_dir = new char[4096];
    output_file = new char[4096];
    input_file = new char[4096];
    affinity_dir = new char[4096];
    
    strcpy(lib_dir, argv[1]);

    rd=atoi(argv[2]);
    nh=atoi(argv[3]);
    Np=atoi(argv[4]);
    
    inc=5;
  
    alleles=new char*[nh];
    for(int i=0; i<nh; i++){
	alleles[i]=new char[4096];
	strcpy(alleles[i], argv[inc+i]);
    }

    alleles_map=new char*[nh];
    for(int i=0; i<nh; i++){
	alleles_map[i]=new char[4096];
	strcpy(alleles_map[i], argv[inc+nh+i]);
    }
    
    strcpy(output_file, argv[inc+2*nh]);
    strcpy(affinity_dir, argv[inc+2*nh+1]);
    strcpy(input_file, argv[inc+2*nh+2]);
  
    ///////////////////////////////
    //Load the output of MixMHCpred
    ///////////////////////////////
    
    ifstream file;
    string line;
    char filename[4096];
    char cline[4096];
    char *pch;
    string p;
    int **peptide; peptide=new int*[Np];
    int *lpep; lpep=new int[Np];
    double **rank; rank=new double*[Np];
    for(int n=0; n<Np; n++){ rank[n]=new double[nh];}
    int l;
    


    sprintf(filename, "%s/../temp/MixMHCpred_%d.txt", lib_dir, rd);
    
    file.open(filename);
    for(int i=0; i<12; i++){
	getline(file, line);
    }
    ct=0;
    
    while (!file.eof()) {
	getline(file, line);
	if(line != ""){
	    strcpy(cline, line.c_str());
	    pch = std::strtok(cline, "\t");
	    l=strlen(pch);
	    peptide[ct]=new int[l];
	    for(int c=0; c<l; c++){
		peptide[ct][c]=map(pch[c]);
	    }
	    lpep[ct]=l;
	    for(int i=0; i<3; i++){
		pch = std::strtok(NULL, "\t");
	    }
	    for(int n=0; n<nh; n++){
		pch = std::strtok(NULL, "\t");
		pch = std::strtok(NULL, "\t");
		rank[ct][n]=atof(pch);
	    }
	    ct++;
	}
	
    }
    file.close();
	
   
    
    /////////////////////////////////////////////
    //Load the coefficients of the neural network
    /////////////////////////////////////////////
    
    int Nnode=5;
    int Nfeature=N+1+Lmax-Lmin+1;
    double *bias; bias=new double[Nnode+1];
    double **weight1; weight1=new double*[Nfeature]; for(int j=0; j<Nfeature; j++){weight1[j]=new double[Nnode];}
    double *weight2; weight2=new double[Nfeature];
    
    sprintf(filename, "%s/bias.txt", lib_dir);
    file.open(filename);
    for(int i=0; i<Nnode+1; i++){
	getline(file, line);
	strcpy(cline, line.c_str());
	bias[i]=atof(cline);
    }
    file.close();

    sprintf(filename, "%s/weight.txt", lib_dir);
    file.open(filename);
    for(int j=0; j<Nfeature; j++){
	getline(file, line);
	strcpy(cline, line.c_str());
	pch = std::strtok(cline, "\t");
	weight1[j][0]=atof(pch);
	for(int i=1; i<Nnode; i++){
	    pch = std::strtok(NULL, "\t");
	    weight1[j][i]=atof(pch);
	}
    }

    for(int i=0; i<Nnode; i++){
	getline(file, line);
	strcpy(cline, line.c_str());
	pch = std::strtok(cline, "\t");
	for(int j=1; j<Nnode; j++){
	    pch = std::strtok(NULL, "\t");
	}
	pch = std::strtok(NULL, "\t");
	weight2[i]=atof(pch);
    }
    file.close();
    
    //////////////////////////////////
    //Load the threshold for the ranks
    //////////////////////////////////

    int Nbin=0;
    
    double *bin; bin=new double[1000];
    int Nrange=5;
    double *range; range=new double[Nrange];
    range[0]=0.001;range[1]=0.01;range[2]=0.1;range[3]=1; range[4]=10;
    for(int j=0; j<Nrange; j++){
	for(int i=1; i<10; i++){
	    bin[Nbin]=log(1.0*range[j]*i);
	    Nbin++;
	}
    }
    bin[Nbin]=log(100);
    Nbin++;
    
    double **rank_thr; rank_thr=new double*[Nbin];
    for(int n=0; n<Nbin; n++) {
	rank_thr[n]=new double[nh];
    }
    
    for(int h=0; h<nh; h++){
	sprintf(filename, "%s/PerRank/%s.txt", lib_dir, alleles_map[h]);
	file.open(filename);
	if(file){
	    for(int n=0; n<Nbin; n++) {
		getline(file, line);
		strcpy(cline, line.c_str());
		pch = std::strtok(cline, "\t");
		pch = std::strtok(NULL, "\t");
		rank_thr[n][h]=atof(pch);
	    }
	    file.close();
	} else {  // If the precomputed ranks are not available for this allele, use those from A0201 - but this should never happen.
	    cout<<"WARNING: missing file of pre-computed ranks for "<<alleles_map[h]<<", using values from A0201"<<endl;
	    sprintf(filename, "%s/PerRank/A0201.txt", lib_dir);
	    file.open(filename);
	    if(file){
		for(int n=0; n<Nbin; n++) {
		    getline(file, line);
		    strcpy(cline, line.c_str());
		    pch = std::strtok(cline, "\t");
		    pch = std::strtok(NULL, "\t");
		    rank_thr[n][h]=atof(pch);
		}
		file.close();
	    } else {
		//This is only useful to compute the threshold for the rank, if PerRank/A0201.txt file does not exist.
		cout<<"Missing file PerRank/A0201.txt"<<endl;
		for(int n=0; n<Nbin; n++) {  
		    rank_thr[n][h]=0.1*(Nbin-n);
		}
	    }
	}
    }
    
    //////////////////////////
    //Compute the PRIME scores
    //////////////////////////

    
    //Precompute the position for each allele and each length
    int ***pos; pos=new int**[nh];
    int **npos; npos=new int*[nh];
    for(int n=0; n<nh; n++){
	pos[n]=new int*[Lmax+1];
	npos[n]=new int[Lmax+1];
	for(int l=Lmin; l<=Lmax; l++){
	    pos[n][l]=new int[Lmax];
	    find_pos(pos[n][l], alleles_map[n], l);
	    ct=0;
	    while(pos[n][l][ct] != 0){
		ct++;
	    }
	    npos[n][l]=ct;
	}
    }

   
    FILE *pFile;
    pFile=fopen(output_file,"w");
    
    
    fprintf (pFile, "####################\n");
    fprintf (pFile, "# Output from PRIME (v2.0)\n");
    fprintf (pFile, "# Alleles: %s",alleles[0]); for(int h=1; h<nh; h++){fprintf (pFile, ", %s", alleles[h]);} fprintf (pFile, "\n");
    fprintf (pFile, "# Affinity predictions: MixMHCpred (%s)\n", affinity_dir);
    fprintf (pFile, "# Input file: %s\n", input_file);
    fprintf (pFile, "#\n");
    fprintf (pFile, "# PRIME is freely available for academic users.\n");
    fprintf (pFile, "# Private companies should contact eauffarth@licr.org at the Ludwig Institute for Cancer Research Ltd for commercial licenses.\n");
    fprintf (pFile, "# To cite PRIME2.0, please refer to:\n");
    fprintf (pFile, "# Gfeller et al. Improved predictions of immunogenicity reveal SARS-Cov-2 CD8 T-cell epitopes, BioRxiv (2022).\n");
    fprintf (pFile, "####################\n");

    fprintf (pFile, "Peptide\t\%%Rank_bestAllele\tScore_bestAllele\t\%%RankBinding_bestAllele\tBestAllele");
    for(int n=0; n<nh; n++){
	fprintf (pFile, "\t\%%Rank_%s\tScore_%s\t\%%RankBinding_%s", alleles[n], alleles[n], alleles[n]);
    }
    fprintf (pFile, "\n");
    
    double fr[N];
    int lp, k;
    int mp;
    double lr,t;
    double *score; score=new double[nh];
    double *rank_PRIME; rank_PRIME=new double[nh];

    int min_pos;
    double min_rank;
    
    for(int i=0; i<Np; i++){
	l=lpep[i];
	for(int n=0; n<nh; n++){

	    //Determine the frequency, includin BLOSUM corrections
	    for(int p=0; p<N; p++){
		fr[p]=0;
	    }
  
	    for(int p=0; p<npos[n][l]; p++){
		mp=pos[n][l][p];
		fr[peptide[i][mp]]=fr[peptide[i][mp]]+1;
	    }
	    for(int j=0; j<N; j++){
		fr[j]=fr[j]/npos[n][l];
	    }

	    //Compute the scores
	    score[n]=bias[Nnode];
	    lr=-log(rank[i][n]);
	    for(int j=0; j<Nnode; j++){
		t=bias[j];
		for(int ii=0; ii<N; ii++){
		    t=t+weight1[ii][j]*fr[ii];
		}
		t=t+weight1[N][j]*lr;  // Encoding the affinity
		t=t+weight1[N+l-Lmin+1][j];  // Encoding the length (one-hot encoding)
		t=1/(1+exp(-t));
		score[n]=score[n]+weight2[j]*t;
	    }
	    score[n]=1/(1+exp(-score[n]));

	    
	    //Compute the %rank
	    
	    k=0;
	    while (score[n]<rank_thr[k][n]) {
		k++;
		if(k==Nbin){
		    break;
		}
	    }
	    if (k==Nbin) {
		rank_PRIME[n]=100;
	    } else if (k==0){
		rank_PRIME[n]=exp(bin[0]);
	    } else {
		//Compute the regression between the two bins
		rank_PRIME[n]=bin[k-1]+(bin[k]-bin[k-1])/(rank_thr[k][n]-rank_thr[k-1][n])*(score[n]-rank_thr[k-1][n]);
		rank_PRIME[n]=exp(rank_PRIME[n]);
	    }
	}

	//Find the minimal value for the ranks across the alleles
	min_rank=1000;
	min_pos=-1;
	for(int n=0; n<nh; n++){
	    if(rank_PRIME[n]<min_rank){
		min_rank=rank_PRIME[n];
		min_pos=n;
	    }
	}
	for(int j=0; j<l; j++){
	    fprintf(pFile, "%c", letter[peptide[i][j]]);
	}
	fprintf(pFile, "\t%.3f\t%.6f\t%.3f\t%s", rank_PRIME[min_pos], score[min_pos], rank[i][min_pos], alleles[min_pos]);
	for (int n=0; n<nh; n++) {
	    fprintf(pFile, "\t%.3f\t%.6f\t%.3f", rank_PRIME[n], score[n], rank[i][n]);
	}
	fprintf(pFile, "\n");
    }
    fclose (pFile);

    
}


void find_pos(int *pos, char *h, int l){

    int ct;

    //Ideally, read a file with these information
    
    std::string s(h);

    ct=0;
    if(s=="B0801" || s=="B1401" || s=="B1402" || s=="B3701" || s=="A6802" || s == "B4701" || s=="H-2-Db"){  // 6
	pos[ct]=3;
	ct++;
	for(int j=5; j<l-1; j++){
	    pos[ct]=j;
	    ct++;
	}
    } else if(s=="A0201" || s=="A0202" || s=="A0207" || s=="A0211"){ 
	for(int j=4; j<l-1; j++){
	    pos[ct]=j;
	    ct++;
	}
    } else if(s=="A0203" || s=="A0204" || s=="A0205" || s=="A0206" || s=="A0220"){ 
	for(int j=4; j<l-4; j++){
	    pos[ct]=4;
	    ct++;
	}
	pos[ct]=l-3;
	ct++;
	pos[ct]=l-2;
	ct++;
    } else if(s=="A2501" || s=="A2601"){ 
	for(int j=4; j<l-4; j++){
	    pos[ct]=4;
	    ct++;
	}
	pos[ct]=l-3;
	ct++;
    } else if (s == "A2902" || s == "B1517" || s == "B5802" || s == "C0602" || s == "C0704"  || s == "C1502" || s == "C1505" || s == "C1602" || s == "C1701" || s == "G0101" || s == "G0103" || s == "G0104"){
	for(int j=4; j<l-3; j++){
	    pos[ct]=j;
	    ct++;
	}
	pos[ct]=l-2;
	ct++;
    } else if(s == "A3201" || s == "B1803" || s == "B3901" || s == "B3905" || s == "B3906" || s == "B3924" || s == "B4601" || s == "B4801" || s == "C0102"){ 
	for(int j=3; j<l-2; j++){
	    pos[ct]=j;
	    ct++;
	}
    } else if(s=="H-2-Kb"){ 
	pos[ct]=3;
	ct++;
	for(int j=6; j<l-1; j++){
	    pos[ct]=j;
	    ct++;
	}
    } else{
	for(int j=3; j<l-1; j++){ 
	    pos[ct]=j;
	    ct++;
	}
    }
}

int map(char c){

    //Should do a binary search, but works fine like this
    
    int p;

    if(c < 'M'){
        if(c=='A'){
	    p=0;
	}else if (c=='C'){
	    p=1;
	}else if (c=='D'){
	    p=2;
	}else if (c=='E'){
	    p=3;
	}else if (c=='F'){
	    p=4;
	}else if (c=='G'){
	    p=5;
	}else if (c=='H'){
	    p=6;
	}else if (c=='I'){
	    p=7;
	}else if (c=='K'){
	    p=8;
	}else if (c=='L'){
	    p=9;
	}
    } else if(c > 'M'){
	if (c=='N'){
	    p=11;
	}else if (c=='P'){
	    p=12;
	}else if (c=='Q'){
	    p=13;
	}else if (c=='R'){
	    p=14;
	}else if (c=='S'){
	    p=15;
	}else if (c=='T'){
	    p=16;
	}else if (c=='V'){
	    p=17;
	}else if (c=='W'){
	    p=18;
	}else if (c=='Y'){
	    p=19;
	}
    } else if(c=='M'){
	p=10;
    }
    return(p);
}
