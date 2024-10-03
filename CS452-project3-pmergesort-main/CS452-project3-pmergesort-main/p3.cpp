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
void pmerge(int *a, int *b, int lasta, int lastb , int * output);

int my_rank;			// my CPU number for this process
int p;					// number of CPUs that we have
MPI_Status status;		// return status for receive
int main (int argc, char * argv[]) {
	
    int k =12;
  	int source;				// rank of the sender
  	int dest;				// rank of destination
  	int tag = 0;			// message number
  	char message[100];		// message itself
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
		cin >> size;
	}
	
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//test n
	//cout <<my_rank<< ", "<< size << endl;
	
	int* myarray = new int[size]; 
	for(int i=0; i<size; i++)
		myarray[i]=0;
	
	//initialize array in p0 and bcast to other processes
	int seed = 9;
    //int seed = 2;
    //int seed = 7;
    srand(seed);
	if(my_rank == 0) {
		for(int j=0; j<size; j++) 
			myarray[j] = rand() % 50;
	cout << "INITIAL ARRAY: " << endl;
	for(int i=0; i<size; i++){
		cout << myarray[i]<< " ";
	}
	cout << endl;
	}
	
	
	MPI_Bcast(&myarray[0], size, MPI_INT, 0, MPI_COMM_WORLD);
	
	//test mergeSort
	mergeSort(myarray, 0, size);
	
	
	//Sorted array
	if(my_rank == 0) {
		cout << endl << endl;
		cout << "Sorted array: " << endl;
		for(int i=0; i<size; i++) 
			cout << myarray[i] << " " ;
		cout << endl;
	}
	
	//shutdown procedures
	delete[] myarray;
	//delete[] arrA;
	//delete[] arrB;
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
        output[indexMerge++] = b[index2++];
    }
    
}

void mergeSort(int * a,int first,int last){
    int m = (first+last)/2;
    int * b = &a[m];
	
	
	if(m==1) {
		if(a[0] > b[0]) {
			int temp = a[1];
			a[1] = a[0];
			a[0] = temp;
		}
		//if(my_rank == 0) {
		//cout << "A: " << a[0] << " " << a[1] << endl;
		//cout << "B: " << b[0] << endl;
		//}
		return;
	}
	
    mergeSort(a,first,m);
    mergeSort(b,first,last-m);
    int* outputarr = new int[last]; 

    pmerge(a, b, m, last-m,outputarr);
    //smerge(a, b, m, last-m,outputarr);
	if(my_rank == 0) {
		//cout << m <<  ", OUTPUT ARRAY: " << endl;
		for(int i = 0; i< last; i++){
			a[i] = outputarr[i];
			//cout << a[i] << " " ;
		}
		//cout << endl;
	}
	MPI_Bcast(&a[0], last, MPI_INT, 0, MPI_COMM_WORLD);
		
}

void pmerge(int *a, int * b, int lasta, int lastb , int * output){
	//STAGE 1 - PARTITIONING
	int selectSize = ceil(((double)lasta)/((double)(log2(lasta))));
	//cout << my_rank << ", selectsize: " << selectSize <<endl;
	int shapesPerP = ceil((double)(selectSize*2)/(double)p); //EX: 10/6
	int shapes = selectSize*2;
	//int shapesPerP = floor((double)(selectSize*2)/(double)p); 
	
	
	int tempp = p;
	if((floor((double)(selectSize*2)/(double)p)) < 1){
		shapesPerP = 1;
		tempp = shapes;
	}
	
	int A_Select[selectSize];
	int B_Select[selectSize];
	
	int aRank[selectSize];
	int bRank[selectSize];
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
	}
	
	
	MPI_Bcast(&A_Select[0], selectSize, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&B_Select[0], selectSize, MPI_INT, 0, MPI_COMM_WORLD);
	
	{
		for(int i=0; i<selectSize; i++){
			int localA_Select = A_Select[i];
			int localB_Select = B_Select[i];
		
			localA_Select=rank1(b, 0, lastb, localA_Select);
			localB_Select=rank1(a, 0, lasta, localB_Select);
			
			A_Select[i] = localA_Select;
			B_Select[i] = localB_Select;
			
		}
	}
	
	int l_bound=0;
	int u_bound=0;
	int i=0;
	int j=0;
	int j1=0;
	int k=0;
	int aSizes[selectSize*2];
	int al_bounds[selectSize*2];
	int au_bounds[selectSize*2];
	int bSizes[selectSize*2];
	int bl_bounds[selectSize*2];
	int bu_bounds[selectSize*2];

	//shapes ---------------------------------------------------------------------------------------
	
	if(my_rank==0){
		
		if(B_Select[i] == 0) {
			while(B_Select[i]==0) {
				aSizes[k] = 0;
				al_bounds[k] = 0;
				au_bounds[k] = 0;
				k++;
				i++;
			}
		}
		j+=log2(lasta);
		j1++;
		
		while(i < selectSize && j1 < selectSize) {
			if(B_Select[i] <= j) {
				u_bound = B_Select[i];
				i++;
			}
		
			else if(B_Select[i] > j) {
				u_bound = j;
				j+=log2(lasta);
				j1++;
			}
			
			
			aSizes[k] = u_bound - l_bound;
			al_bounds[k] = l_bound;
			au_bounds[k] = u_bound;
			
			k++;
			
			//set new l_bound
			l_bound = u_bound; 
		}
		
		while(i < selectSize) {
			u_bound = B_Select[i];
			i++;
			
			aSizes[k] = u_bound - l_bound;
			al_bounds[k] = l_bound;
			au_bounds[k] = u_bound;	
			
			k++;
			
			//set new l_bound
			l_bound = u_bound;
		}
		
		while(j1 < selectSize) {
			
			u_bound = j;
			j+=log2(lasta);
			j1++;
			
			aSizes[k] = u_bound - l_bound;
			al_bounds[k] = l_bound;
			au_bounds[k] = u_bound;	
			
			k++;
			//set new l_bound
			l_bound = u_bound;
		}
		
		u_bound = lasta;
		
		aSizes[k] = u_bound - l_bound;
		al_bounds[k] = l_bound;
		au_bounds[k] = u_bound;
	} //end if (my_rank == 0)
	
	//receive bounds data
	MPI_Bcast(&aSizes[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&au_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&al_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
	

	//Calculate for b
	
	if(my_rank==0){
		i=j=k=j1=l_bound=u_bound=0;
		
		if(A_Select[i] == 0) {
			while(A_Select[i]==0) {
				bSizes[k] = u_bound - l_bound;
				bl_bounds[k] = l_bound;
				bu_bounds[k] = u_bound;
				k++;
				i++;
			}
			
		}
		j+=log2(lasta);
		j1++;
		
		while(i < selectSize && j1 < selectSize) {
		
			if(A_Select[i] <= j) {
				u_bound = A_Select[i];
				i++;
			}
		
			else if(A_Select[i] > j) {
				u_bound = j;
				j+=log2(lasta);
				j1++;
			}
			
			
			bSizes[k] = u_bound - l_bound;
			bl_bounds[k] = l_bound;
			bu_bounds[k] = u_bound;
			
			k++;
			//set new l_bound
			l_bound = u_bound; 
		
		}
		
		while(i < selectSize) {
			u_bound = A_Select[i];
			i++;
			
			bSizes[k] = u_bound - l_bound;
			bl_bounds[k] = l_bound;
			bu_bounds[k] = u_bound;	
			
			k++;
			//set new l_bound
			l_bound = u_bound;
		}
		
		while(j1 < selectSize) {
			
			u_bound = j;
			j+=log2(lasta);
			j1++;
			
			bSizes[k] = u_bound - l_bound;
			bl_bounds[k] = l_bound;
			bu_bounds[k] = u_bound;	
			
			k++;
			//set new l_bound
			l_bound = u_bound;
		}
		
		u_bound = lasta;
		
		bSizes[k] = u_bound - l_bound;
		bl_bounds[k] = l_bound;
		bu_bounds[k] = u_bound;
		
	} //end if 0
	

	//receive bounds data
	MPI_Bcast(&bSizes[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&bu_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&bl_bounds[0], selectSize*2, MPI_INT, 0, MPI_COMM_WORLD);
	
	/*cout <<my_rank<<", al_bounds: ";
	for(int i=0; i<selectSize*2; i++)
		cout << al_bounds[i] << " " ;
	cout << endl;
	cout << my_rank << ", au_bounds: ";
	for(int i=0; i<selectSize*2; i++)
		cout << au_bounds[i] << " " ;
	cout << endl;
	cout <<my_rank<<", bl_bounds: ";
	for(int i=0; i<selectSize*2; i++)
		cout << bl_bounds[i] << " " ;
	cout << endl;
	cout <<my_rank<<", bu_bounds: ";
	for(int i=0; i<selectSize*2; i++)
		cout << bu_bounds[i] << " " ;
	cout << endl;*/
	
//shapes
	int shapesPerP1 = shapesPerP;
	if(my_rank == (tempp-1)) {
		
		if((shapesPerP*(tempp-1))==shapes){
			shapesPerP1 = -1;
		}
		else if((selectSize*2)%tempp != 0) {
			shapesPerP1 = (selectSize*2)%shapesPerP;
		}
		
	}
	
	if(my_rank < tempp) {
			
		for(int i=0; i<=shapesPerP1; i++) {
			int starta = al_bounds[(my_rank*shapesPerP)+i];
			int startb = bl_bounds[(my_rank*shapesPerP)+i];
			int lasta1 = au_bounds[(my_rank*shapesPerP)+i];
			int lastb1 = bu_bounds[(my_rank*shapesPerP)+i];
			int sizea = aSizes[(my_rank*shapesPerP)+i];
			int sizeb = bSizes[(my_rank*shapesPerP)+i];
			int size = sizea + sizeb;
		
			smerge(&a[starta], &b[startb], sizea, sizeb, &output[0]);
				
			MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&output[0], size, MPI_INT, 0, 0, MPI_COMM_WORLD);
			
			if(i==shapesPerP1-1){
				i++;
			}
		}
	}
		if(my_rank==0) {
			int placeholder = 0;
			int recvSize = 0;
			for(int i=0; i<tempp; i++) {
				//cout << "I: " << i << endl;
				if(tempp == p) {
					if(i == (tempp-1)) {
						
							if((shapesPerP*(tempp-1))==shapes){
								shapesPerP = -1;
							}
						
							else if((selectSize*2)%tempp != 0) {
								shapesPerP = (selectSize*2)%shapesPerP;
							}
						
					}
				}
				
				for(int j=0; j<=shapesPerP; j++) {
					MPI_Recv(&recvSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(&output[placeholder], recvSize, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					placeholder += recvSize;
					if(j==shapesPerP-1){
						j++;
					}
					
				}
			}
		
		} //end if my_rank==0
	
}//end pmerge