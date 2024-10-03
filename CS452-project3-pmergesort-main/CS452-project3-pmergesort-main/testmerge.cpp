#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include <math.h>
using namespace std;

// New compile and run commands for MPI!
// mpicxx -o blah testmerge.cpp
// mpirun -q -np 4 blah

int rank1(int * a, int first, int last, int valToFind);
void smerge(int *a, int *b, int lasta, int lastb, int * output);
void mergeSort(int * a, int first, int last);
void pmerge(int *a, int *b, int lasta, int lastb , int * output );
int my_rank;			// my CPU number for this process
  	int p;					// number of CPUs that we have

int main (int argc, char * argv[]) {
	
    int k =12;
  	int source;				// rank of the sender
  	int dest;				// rank of destination
  	int tag = 0;			// message number
  	char message[100];		// message itself
  	MPI_Status status;		// return status for receive
	int evenmax;
	int oddmax;
  	// Start MPI
  	MPI_Init(&argc, &argv);

  	// Find out my rank!
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  	// Find out the number of processes!
  	MPI_Comm_size(MPI_COMM_WORLD, &p);


	int size = 0;
	
	if(my_rank == 0) {
		cout<<"how big is the array: " << endl;
		//change later 
		//size = 64;
		cin >> size;
	}
	
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//test n
	cout << size << endl;
	
	int* myarray = new int[size]; 
	for(int i=0; i<size; i++)
		myarray[i]=0;
	
	//initialize array in p0 and bcast to other processes
	//int seed = 9;
    //int seed2 = 2;
    int seed3 = 7;
    srand(seed3);
	if(my_rank == 0) {
		for(int j=0; j<size; j++)
			myarray[j] = rand() % 50;
	}
	
	MPI_Bcast(&myarray[0], size, MPI_INT, 0, MPI_COMM_WORLD);

	
	//test bcast
	/*for(int i=0; i<size; i++){
		cout << myarray[i]<< " ";
	}*/
	cout << endl;
	
	int * arrA = new int[size/2];
	int * arrB = new int[size/2];
	
	//test pmerge sranks
	if(my_rank == 0) {
		
		
		for(int i=0; i<size/2; i++)
			arrA[i] = myarray[i];
		
		
		for(int i=size/2; i<size; i++)
			arrB[i-(size/2)] = myarray[i];
		
		mergeSort(arrA, 0, size/2);
		mergeSort(arrB, 0, size/2);
		
		/*for(int i=0; i<size/2; i++){
			cout << "A: " << arrA[i]<< " ";
		}*/
		
		for(int i=0; i<size/2; i++){
			cout << "B: " << arrB[i]<< " ";
		}
		cout << endl << endl;
		
		
		
		//int rankval = rank1(myarray, 0, size, 50);
		//cout << "RANK VALUE of 50: " << rankval << endl;
	}
	
	MPI_Bcast(&arrA[0], size/2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&arrB[0], size/2, MPI_INT, 0, MPI_COMM_WORLD);
	
	pmerge(arrA, arrB, size/2, size/2,myarray); 

	//cout << "RANK VALUE OF 24: " << rank1(arrB, 0, size/2, 24) << endl;
	
	//shutdown procedures
	delete[] myarray;
	delete[] arrA;
	delete[] arrB;
	MPI_Finalize();
	return 0;   
}


int rank1(int * a, int first, int last, int valToFind){
	
	//create aliases for parts of arrays
	int firstpartStart = a[first];
	int secondhalfStart = a[last/2];
	
	//find number of items less than valToFind
	int n=0;
	
	if(valToFind >= secondhalfStart)
		n = last/2;
	
	while((valToFind > a[n]) && (n < last)) {
		n++;
	}
	
	return n;
}


void smerge(int * a, int * b, int lasta, int lastb, int * output){
    int index1 = 0;
    int index2 =0;
    int indexMerge =0;

    //cout << "MERGE lasta " << lasta << endl;
    //cout << "MERGE lastb " << lastb << endl;
	if((lasta+lastb)!=0){
    while (index1 < lasta && index2 < lastb) {
        if (a[index1] <= b[index2]) {
            output[indexMerge++] = a[index1++];              
                   
        }
        else {
            output[indexMerge++] = b[index2++];             
        }
        
    }

    while (index1 < lasta ){
        output[indexMerge++]= a[index1++];
    }
    
   while (index2 < lastb ){
        //cout <<"mergeindex  " <<indexMerge<< endl;
        //cout<<"index2 " << index2 <<endl;
        //cout << "MERGE lastb " << lastb << endl;
        output[indexMerge++] = b[index2++];
    }
	}
    /*for(int i =0; i<indexMerge; i++){
        a[i] = output[i];
    }*/
	
}

void mergeSort(int * a,int first,int last){
    int m = (first+last)/2;
    int * b = &a[m];

	//temp
	/*if(last<= 32){
		return;
	}*/
	
    if(first>=m){
        return;
    }

    //cout << "let it end " <<first <<endl;
    mergeSort(a,first,m);
    mergeSort(b,first,last-m);
    int* outputarr = new int[last]; 

   smerge(a, b, m, last-m,outputarr);
    for(int i = 0; i< last; i++){
        a[i] = outputarr[i];
    }
}

void pmerge(int *a, int * b, int lasta, int lastb , int * output = NULL){
		MPI_Status status;
	//STAGE 1 - PARTITIONING
	int selectSize = ceil(((double)lasta)/((double)(log2(lasta))));
	cout << "selectsize: " << selectSize <<endl;
		int shapesPerP = (selectSize*2)/p;
		cout <<"shapesperP: "<<shapesPerP<<endl;
	//SRANKA = A(1), A(1 + log n), A(1 + 2 log n). . . A(n + 1 − log n) (change to 0 indexing)
	//a and b same n so use same selectSize for both
	
	
	//use selectSize number of ps in this
		
	int * A_Select = new int[selectSize];
	int * B_Select = new int[selectSize];
	
	int * aRank = new int[selectSize];
	int * bRank = new int[selectSize];
	if(my_rank == 0) {
		//filling A_Select and B_Select
		int j=0;
		for(int i=0; i<selectSize; i++) {
			A_Select[i] = a[j];
			B_Select[i] = b[j];
			aRank[i] = j;
			bRank[i] = j;
			j += ceil((double)(log2(lasta)));
			
		
		}
		
	/*
		cout << "INITIAL SELECTA/B" << endl;
		for(int i=0; i<selectSize; i++)
			cout << A_Select[i] << " ";
		cout << endl;
			
		for(int i=0; i<selectSize; i++)
			cout << B_Select[i] << " ";
		cout << endl;
	*/
	}
	
	
	//rank in parallel
	
	MPI_Bcast(&A_Select[0], selectSize, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&B_Select[0], selectSize, MPI_INT, 0, MPI_COMM_WORLD);
	
	//THIS IS SRANK A AND B MOTHERFUCKERS
	{
		for(int i=0; i<selectSize; i++){
		int localA_Select = A_Select[i];
		int localB_Select = B_Select[i];
		
		//cout << "MYRANK: " << my_rank << "LOCALA: " << localA_Select << endl;
		localA_Select=rank1(b, 0, lastb, localA_Select);
		//cout << "MYRANK: " << my_rank << "LOCALB: " << localB_Select << endl;
		localB_Select=rank1(a, 0, lasta, localB_Select);
		if(my_rank==0){
		cout << "MYRANK: " << my_rank << "sendLOCALB: " << localB_Select << endl;
		cout << "MYRANK: " << my_rank << "sendLOCALA: " << localA_Select << endl;
		}
			A_Select[i] = localA_Select;
			B_Select[i] = localB_Select;
			
		}
		
	}
	//cout << "MYRANK: " << my_rank << endl;
	
	
	//STAGE 1B - put srank/select vals into WIN

	int WIN[lasta+lastb];
	
	
	for(int i=0; i<lasta+lastb; i++)
		WIN[i]=0;
	
	{
		int selectvalA=0;
		int selectvalB=0;
		int j=0;
		for(int i=0; i<selectSize; i++){

			selectvalA = A_Select[i];
			selectvalB = B_Select[i];
			WIN[selectvalA+j] = a[j];
			WIN[selectvalB+j] = b[j];
			j += ceil((double)(log2(lasta)));
		}
}
//ranka = 1,4,5 a's position in a
//srankb = 0,2,8 b's rank in a
//lb=0 ub=1 lb=2 ub=7 lb=8 ub=
	int l_bound=0;
	int u_bound=0;
	int srankbcounter=0;
	int srankacounter=0;
	int j=0;
	int arankcounter=0;
	int brankcounter=0;
	int k=0;
	int aSizes[selectSize*2];
	int al_bounds[selectSize*2];
	int au_bounds[selectSize*2];
	int bSizes[selectSize*2];
	int bl_bounds[selectSize*2];
	int bu_bounds[selectSize*2];

	//shapes ---------------------------------------------------------------------------------------
	if(my_rank==0){
	
	cout << "A_Select: " << endl;
	for(int i=0; i<selectSize; i++)
		cout << A_Select[i] << " ";
	cout << endl;
	cout << "B_Select: " << endl;
	for(int i=0; i<selectSize; i++)
		cout << B_Select[i] << " ";
	cout << endl;

		cout <<"arank: ";
		for(int i=0; i<selectSize; i++)
		cout << aRank[i] << " ";
	cout << endl;

		cout <<"bRank: " <<endl;
		for(int i=0; i<selectSize; i++)
		cout << bRank[i] << " ";
	cout << endl;

	}//subarray of albounds between [(n-1)…(n)] of all the contents, move to win by upper and lower bounds
	//B_Select is srankb
	//A_select is sranka
	{
		while(srankbcounter < selectSize && arankcounter < selectSize) {//all of a's lower bounds

			if(B_Select[srankbcounter] <aRank[arankcounter]) {
			al_bounds[l_bound] = B_Select[srankbcounter];
			/*if(my_rank==0){
			cout<< "albounds in first if " << al_bounds[l_bound] <<endl;
			cout <<"srankb counter " << srankbcounter <<endl;
			


			
			cout << my_rank << " rank is srankb lower bound: "<<l_bound<<endl;
			}*/
			srankbcounter++;
			l_bound++;
			}
			else {
				al_bounds[l_bound]=aRank[arankcounter];
				
				/*if(my_rank==0){
					cout<< "albounds in second  " << al_bounds[l_bound] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<l_bound<<endl;
				cout <<"arank counter " << arankcounter <<endl;


				}*/
				arankcounter++;
				l_bound++;
				}
			}
			{
			while(arankcounter < selectSize){
				al_bounds[l_bound]=aRank[arankcounter];
				

				/*if(my_rank==0){
					cout<< "albounds in third  " << al_bounds[l_bound] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<l_bound<<endl;
				cout <<"arank counter " << arankcounter <<endl;


				}*/
				arankcounter++;
				l_bound++;
				}

			}
			{
			while(srankbcounter < selectSize){
				al_bounds[l_bound]=B_Select[srankbcounter];
			/*if(my_rank==0){
				cout<< "albounds in fourth  " << al_bounds[l_bound] <<endl;
				cout << my_rank << " rank is srankb lower bound: "<<l_bound<<endl;
				cout <<"srankb counter " << srankbcounter <<endl;

				}*/


				srankbcounter++;
				l_bound++;
				}

			}
			if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
				cout  <<" A lowerbounds: ";
			for(int i=0; i<selectSize*2; i++){
				cout<< al_bounds[i] << " ";
			}
			cout<<endl;
			MPI_Bcast(&al_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
			}

	}
//b lower bound ------------------------------------------------
	{
		l_bound=0;
		while(srankacounter < selectSize && brankcounter < selectSize) {//all of a's lower bounds
			if(A_Select[srankacounter] <bRank[brankcounter]) {
				bl_bounds[l_bound] = A_Select[srankacounter];
				/*if(my_rank==0){
				cout<< "albounds in first if " << al_bounds[l_bound] <<endl;
				cout <<"srankb counter " << srankbcounter <<endl;
				cout << my_rank << " rank is srankb lower bound: "<<l_bound<<endl;
				}*/
				srankacounter++;
				l_bound++;
			}
			else {
				bl_bounds[l_bound]=bRank[brankcounter];
				/*if(my_rank==0){
					cout<< "albounds in second  " << al_bounds[l_bound] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<l_bound<<endl;
				cout <<"arank counter " << arankcounter <<endl;
				}*/
				brankcounter++;
				l_bound++;
			}
		}
		
		while(brankcounter < selectSize){
			bl_bounds[l_bound]=bRank[brankcounter];
				/*if(my_rank==0){
					cout<< "albounds in third  " << al_bounds[l_bound] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<l_bound<<endl;
				cout <<"arank counter " << arankcounter <<endl;


				}*/
				brankcounter++;
				l_bound++;
		}	
		while(srankacounter < selectSize){
			bl_bounds[l_bound]=A_Select[srankacounter];
			/*if(my_rank==0){
				cout<< "albounds in fourth  " << al_bounds[l_bound] <<endl;
				cout << my_rank << " rank is srankb lower bound: "<<l_bound<<endl;
				cout <<"srankb counter " << srankbcounter <<endl;

			}*/
			srankacounter++;
			l_bound++;
		}

			
		if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
			cout  <<" B lowerbounds: ";
			for(int i=0; i<selectSize*2; i++){
				cout<< bl_bounds[i] << " ";
			}
			cout<<endl;
			MPI_Bcast(&bl_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	//Upper bounds -- A 
	{
		int indexAlBounds = 0;
		srankbcounter= 0;
		arankcounter = 0;
		while(indexAlBounds != selectSize*2){
			if(al_bounds[indexAlBounds] ==0){ // do this for if the lower bounds ==0 then upper bound = 0 -- i might be wrong for this
				au_bounds[u_bound] = 0;
				u_bound++;
				indexAlBounds++;
				//have to increment corrisponding arankbcounter and arank counter 
				if(aRank[arankcounter] == 0){
					arankcounter++;
				}
				else{
					srankbcounter++;
				}
			}
			else{ //do this if alBounds[i] != 0
				if(al_bounds[indexAlBounds] == B_Select[srankbcounter]){
					au_bounds[u_bound] = al_bounds[indexAlBounds-1]+1;
					u_bound++;
					indexAlBounds++;
					srankbcounter++;
				}
				else{// if(al_bounds[indexAlBounds] == aRank[arankcounter]){
					au_bounds[u_bound] = al_bounds[indexAlBounds-1]+1;
					u_bound++;
					indexAlBounds++;
					arankcounter++;
				}
				
			}
			
		}
		if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
			cout  <<" A upperbounds: ";
			for(int i=0; i<selectSize*2; i++){
				cout<< au_bounds[i] << " ";
			}
			cout<<endl;
			MPI_Bcast(&au_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
		}
		

	}
	//Upper bounds - B
	{
		int indexBlBounds = 0;
		u_bound =0;
		srankacounter = 0;
		brankcounter = 0;
		while(indexBlBounds != selectSize*2){
			if(bl_bounds[indexBlBounds] ==0 ){ // do this for if the lower bounds ==0 then upper bound = 0 -- i might be wrong for this
				//have to increment corrisponding arankbcounter and arank counter 
				if(bRank[brankcounter] == 0){
					bu_bounds[u_bound] = 0;
					brankcounter++;
				}
				else if(A_Select[arankcounter] ==0){
					bu_bounds[u_bound] = 0;
					srankacounter++;
				}
				u_bound++;
				indexBlBounds++;
			}
			else{ //do this if alBounds[i] != 0
				if(bl_bounds[indexBlBounds] == A_Select[srankacounter]){
					bu_bounds[u_bound] = bl_bounds[indexBlBounds-1]+1;
					u_bound++;
					indexBlBounds++;
					srankacounter++;
				}
				else{// if(al_bounds[indexAlBounds] == aRank[arankcounter]){
					bu_bounds[u_bound] = bl_bounds[indexBlBounds-1]+1;
					u_bound++;
					indexBlBounds++;
					brankcounter++;
				}
				
			}
			
		}
		if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
			cout  <<" B upperbounds: ";
			for(int i=0; i<selectSize*2; i++){
				cout<< bu_bounds[i] << " ";
			}
			cout<<endl;
			MPI_Bcast(&bu_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
		}
		
		

	}
	//a and b sizes
		/*if (my_rank==0){
			cout <<" asizes: ";
			for(int i=0; i<selectSize*2;i++){
				aSizes[i]= al_bounds[i] - au_bounds[i];
				cout << aSizes[i]<<" ";
			
			}
			cout<<endl;
			MPI_Bcast(&aSizes[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
		}

		
		if (my_rank==0){
			cout <<" bsizes: ";
			for(int i=0; i<selectSize*2;i++){
				bSizes[i]= bl_bounds[i] - bu_bounds[i];
				cout << bSizes[i]<<" ";
			
			}
			cout<<endl;
			MPI_Bcast(&bSizes[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
		}*/

	for(int u = my_rank; u < (2*selectSize); u+=p){




	}


for(int i=0; i<lasta+lastb;i++){
				
				output[i]=0;
			
			}
	//receive and sort
		//1. bcast aSizes
		//2. create local arrays based on aSizes
		//3. receive into local arrays
//	int *  outputarrs = new int[shapesPerP][];
//cout << "merge a: " <<endl;
/*for(int u = my_rank; u < (selectSize); u+=p){
cout << b[bl_bounds[(int)u-1]+1] << " to: " << b[bl_bounds[(int)u]] << " ";

if(a[al_bounds[(int)u+1]] ==0) {

a[al_bounds[(int)u+1]] =1;


}
if(b[bl_bounds[(int)u+1]] ==0) {

b[bl_bounds[(int)u+1]] =1;


}


}*/
int shapesize = 0;

if(my_rank==0){

}
else{
	smerge(&a[al_bounds[my_rank-1]+1], &b[bl_bounds[my_rank-1]+1], al_bounds[my_rank], bl_bounds[my_rank], &WIN[0]);
//first a first b last a lastb output
	shapesize =(al_bounds[my_rank]-al_bounds[my_rank-1]+1)+(bl_bounds[my_rank]-bl_bounds[my_rank-1]+1);
	cout << shapesize << " big rank: " <<my_rank <<endl;
}


//MPI_Allreduce( const void* sendbuf , void* recvbuf , MPI_Count count , MPI_Datatype datatype , MPI_Op op , MPI_Comm comm);

MPI_Allreduce(output, WIN,shapesize , MPI_INT, MPI_SUM, MPI_COMM_WORLD);
cout<<"test output: " << endl;
for(int i=0; i<lasta+lastb;i++){
				
				cout << output[i]<<" ";
			
			}

















/*int shapesPerP1 = shapesPerP;
	if(my_rank == (p-1)) {
		if((selectSize*2)%p != 0) {
			shapesPerP1 = (selectSize*2)%shapesPerP;
		}
	}




		
//shapes
	if(my_rank < p) {		 
		for(int i=0; i<=shapesPerP1; i++) {
			int lasta1 = al_bounds[(my_rank*shapesPerP)+i];
			int lastb1 = bl_bounds[(my_rank*shapesPerP)+i];
			int starta = au_bounds[(my_rank*shapesPerP)+i];
			int startb = bu_bounds[(my_rank*shapesPerP)+i];
			int sizea = aSizes[(my_rank*shapesPerP)+i];
			int sizeb = bSizes[(my_rank*shapesPerP)+i];
			//cout << my_rank << ", OUTPUT - starta: " << starta << " startb: " << startb << " lasta: " << lasta1 << " lastb: " << lastb1 << " sizea: " << sizea << " sizeb: " << sizeb << endl;	
			//cout << my_rank << ", A[STARTA]: " << a[starta] << " B[STARTB]: " << b[startb] << endl;
			int size = sizea + sizeb;
			int output[size];
		
			smerge(&a[starta], &b[startb], sizea, sizeb, &output[0]);
				
			cout << my_rank << ", SORTED OUTPUT: ";
			for(int i=0; i<size; i++)
				cout << output[i] << " ";
			cout << endl;
			
			MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&output[0], size, MPI_INT, 0, 0, MPI_COMM_WORLD);
			
			if(i==shapesPerP1-1){
				i++;
			}
		}
		//cout << my_rank << ", made it here (487)" << endl;
		//end if my_rank==0
	}*/
		//cout << my_rank << ", SORTED OUTPUT: ";
		//for(int i=0; i<sizea+sizeb; i++)
		//	cout << output[i] << " ";
		//cout << endl;
		//for loop add 1 multiply shapes per p
		
	



		
		/*if(my_rank == 0) {
			int placeholder = 0;
			for(int i=0; i<p; i++) {
				size = aSizes[my_rank*2] + bSizes[my_rank*2];
				size2 = aSizes[(my_rank*2)+1] + bSizes[(my_rank*2)+1];
				MPI_Recv(&WIN[placeholder], size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
				cout << "Placeholder: " << placeholder << endl;
				placeholder += size;
				cout << "Placeholder: " << placeholder << endl;
				MPI_Recv(&WIN[placeholder], size2, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
				placeholder += size2;
			
			}
			cout << "SORTED WIN ARRAY" << endl;
			cout << lasta + lastb << endl;
			for(int i=0; i<(lasta+lastb); i++)
				cout << WIN[i] << " ";
			cout << endl;
			
				
		}*/
		//cout << my_rank << ", SORTED OUTPUT2: ";
		//for(int i=0; i<sizea+sizeb; i++)
		//	cout << output2[i] << " ";
		//cout << endl;
	
	
		
		
	
	
}

