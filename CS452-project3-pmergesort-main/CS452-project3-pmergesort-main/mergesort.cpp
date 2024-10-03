#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
using namespace std;

// New compile and run commands for MPI!
// mpicxx -o blah file.cpp
// mpirun -q -np 4 blah
int rank1(int * a, int first, int last, int valToFind);
void smerge(int *a, int *b, int lasta, int lastb, int * output = NULL);
void mergeSort(int * a, int first, int last);
void pmerge(int *a, int *b, int lasta, int lastb, int * output = NULL);
void mergeSortp(int *a, int first, int last);


int main (int argc, char * argv[]) {
	int my_rank;			// my CPU number for this process
  	int p;					// number of CPUs that we have
    int k =12;
  	int source;				// rank of the sender
  	int dest;				// rank of destination
  	int tag = 0;			// message number
  	char message[100];		// message itself
  	MPI_Status status;		// return status for receive
	int evenmax;
	int oddmax;
	int num_merge=0; //number of merges 
  	// Start MPI
  	MPI_Init(&argc, &argv);

  	// Find out my rank!
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  	// Find out the number of processes!
  	MPI_Comm_size(MPI_COMM_WORLD, &p);


	int size = 0;
	
	if(my_rank == 0) {
		//cout<<"how big is the array: " << endl;
		//change later 
		size = 64;
		//cin >> size;
	}
	
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//test n
	cout << size << endl;
	
	int* myarray = new int[size]; 
	for(int i=0; i<size; i++)
		myarray[i]=0;
	
	//initialize array in p0 and bcast to other processes
	int seed = 9;
    int seed2 = 2;
    int seed3 = 7;
    srand(seed);
	if(my_rank == 0) {
		for(int j=0; j<size; j++)
			myarray[j] = rand() % 50;
	}
	
	//MPI_Bcast(&myarray[0], size, MPI_INT, 0, MPI_COMM_WORLD);

	

    /*if (rank == 0) {
        // send data array to process 1
        MPI_Send(&myarray, 5, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else if (rank == 1) {
        // receive data array from process 0
        MPI_Recv(&myarray, 5, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // print received data
        for (int i = 0; i < 5; i++) {
            cout << myarray[i] << " ";
        }
        cout << endl;
    }*/
	
	//test bcast
	if(my_rank == 0){
		for(int i=0; i<size; i++){
			cout <<"initial: " << myarray[i]<<  endl;
		}
	}
	
	//cout << endl;
	
	//test rank
	/*
	if(my_rank == 0) {
		
		mergeSort(myarray, 0, size);
		
		for(int i=0; i<size; i++){
			cout << myarray[i]<< " ";
		}
		cout << endl;
		
		int rankval = rank1(myarray, 0, size, 50);
		cout << "RANK VALUE of 50: " << rankval << endl;
	}*/
	//end of smerge
	//gather first then sort 
	// pmerge time
	//size = 64;
	int sRankSize = (size/p);
	int * RankA = new int[size/2];
	int * RankB = new int[size/2];
	int * WIN = new int[size];
	
	int * SRankA = new int[sRankSize];
	int * SRankB = new int[sRankSize];
	int * tempArray = new int[size/p];

	int * outputarr = new int[16];
	int testa[8] = {2,4,6,8,10,12,14,16};
	int testb[8] = {1,3,5,7,9,11,13,15};
	MPI_Scatter(&myarray[0], size/p, MPI_INT, tempArray, size/p, MPI_INT, 0, MPI_COMM_WORLD);

	//testing scatter
	cout << "Current Processor: " << my_rank << ". Temp Array contents: ";
	for (int i = 0; i < size/p; i++){
		cout << tempArray[i] << " " <<endl;
	}
	cout <<endl;

	cout << "Current Processor: " << my_rank <<endl;
	mergeSort(tempArray, 0, size/p);
	if(p%2 ==0){
		if(my_rank%2!= 0){
			MPI_Send(tempArray, size/p, MPI_INT, tempArray, my_rank-1, tag, MPI_COMM_WORLD); //send temp array from odd processor to even (odd p -1)
		}
	}

	
	cout << "is the mergesort failing?" <<endl;

	//localSRankA & localSRankB are the SRankA and SRankB for each Processor
	int *localSRankA = new int[size/p];
	int *localSRankB= new int[size/p];
	
	//if even number of processors
	if(p%2== 0){
		if(my_rank%2 == 0){
			// int * b = &a[m];
			int * btemp = new int[size/p];
			MPI_Recv(btemp, size/p, MPI_INT, my_rank+1, tag, MPI_COMM_WORLD, &status);// recv temp array from odd processor(even+1) in btemp in current P
			//calculate local SRankA
			cout << "Local SRankA: ";
			for(int i =0; i < size/p; i ++){
				localSRankA[i] = rank1(btemp,0, size/p, tempArray[i]);
				cout << localSRankA[i] << " ";
			}
			cout << endl;
			//calculate local SRANKB
			int * atemp = &tempArray[0];
			cout << "Local SRankB: ";
			for(int i = 0; i <size/p; i ++){
				localSRankB[i] = rank1(atemp,0,size/p, btemp[i]);
				cout << localSRankB[i] << " ";
			}
			cout << endl;

			pmerge(localSRankA, localSRankB, size/p, size/p, WIN);

			//storing in Main SRanks;

		}


	}
	
	
	/*
	for(int i = 0;i <= p;i=i+2){
		int selecta = i;
		int selectb = i+1;
		if(my_rank == selecta){
			SRankA = tempArray;
		}
		if(my_rank == selectb){
			SRankB = tempArray;
		}
	}*/
	



	
	//shutdown procedures
	delete[] myarray;
	delete[] RankA;
	delete[] RankB;
	delete[] tempArray;
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
    for(int i =0; i<indexMerge; i++){
        a[i] = output[i];
    }

}

void mergeSort(int * a,int first,int last){
    int m = (first+last)/2;
    int * b = &a[m];

	//temp
	/*if(last<= 32){
		return;
	}
	*/
    if(first>=m){
        return;
    }

    //cout << "let it end " <<first <<endl;
    mergeSort(a,first,m);
    mergeSort(b,first,last-m);
    int* outputarr = new int[last]; 

	//smerge
    smerge(a, b, m, last-m,outputarr);
    
	/*for(int i = 0; i< last; i++){
        a[i] = outputarr[i];
    }*/

	//pmerge 
	//pmerge(a, b, m, last-m, outputarr);
	for (int i = 0; i< last;i++){
		a[i] = outputarr[i];
	}


}
void mergeSortp(int *a, int first, int last){
	int m = (first+last)/2;
	int * b =&a[m];

	if(last>=m){
		return;
	}

	int * outp = new int[last];


}

void pmerge(int *a, int * b, int lasta, int lastb, int * output){
	int outsize = (lasta-1)+(lastb-1);
	for(int i = 0; i < outsize; i++ ){
		if(a[i]<b[i]){
			output[i] = b[i];
		}
	}

}