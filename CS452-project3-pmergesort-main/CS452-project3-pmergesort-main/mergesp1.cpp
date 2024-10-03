#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include "mpi.h"
#include <time.h>
using namespace std;

// New compile and run commands for MPI!
// mpicxx -o blah file.cpp
// mpirun -q -np 4 blah
//control c to quit operation
int p;
int my_rank;
int s;

int Rank(int * a, int first, int last, int valToFind);
void smerge(int *a, int *b, int lasta, int lastb, int * output = NULL);
void mergeSort(int * a, int first, int last);
void pmerge(int *a, int *b, int lasta, int lastb, int * output = NULL);

int main (int argc, char * argv[]) {
	//int my_rank;			// my CPU number for this process
  	//int p;					// number of CPUs that we have
    int k =12;
  	int source;				// rank of the sender
  	int dest;				// rank of destination
  	int tag = 0;			// message number
  	char message[100];		// message itself
  	MPI_Status status;		// return status for receive

  	// Start MPI
  	MPI_Init(&argc, &argv);

  	// Find out my rank!
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  	// Find out the number of processes!
  	MPI_Comm_size(MPI_COMM_WORLD, &p);


	s = 0;
	
	if(my_rank == 0) {
		//cout<<"how big is the array: " << endl;
		//change later 
		s = pow(2,p);
		//cin >> size;
	}
	
	MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//test n
	//cout << s << endl;
	
	int* myarray = new int[s]; 
	for(int i=0; i<s; i++)
		myarray[i]=0;
	
	//initialize array in p0 and bcast to other processes
	int seed = 9;
    int seed2 = 2;
    int seed3 = 7;
    srand(seed);
	if(my_rank == 0) {
		for(int j=0; j<s; j++)
			myarray[j] = rand() % 50;
	}
	
	MPI_Bcast(&myarray[0], s, MPI_INT, 0, MPI_COMM_WORLD);

	//test initial arrau 
    
	if(my_rank == 0){
        cout << "Initial array: ";
		for(int i=0; i<s; i++){
			cout << myarray[i] << " ";
		}
        cout << endl;
	}
int * arrA = new int[s/2];
	int * arrB = new int[s/2];
	
	//test pmerge sranks
	if(my_rank == 0) {
		
		
		for(int i=0; i<s/2; i++)
			arrA[i] = myarray[i];
		
		
		for(int i=s/2; i<s; i++)
			arrB[i-(s/2)] = myarray[i];
		
		//mergeSort(arrA, 0, s/2);
		//mergeSort(arrB, 0, s/2);
		
		/*for(int i=0; i<size/2; i++){
			cout << "A: " << arrA[i]<< " ";
		}*/
		
		for(int i=0; i<s/2; i++){
			cout << "B: " << arrB[i]<< " ";
		}
		cout << endl << endl;
		
		
		
		//int rankval = rank1(myarray, 0, size, 50);
		//cout << "RANK VALUE of 50: " << rankval << endl;
	}
	
	MPI_Bcast(&arrA[0], s/2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&arrB[0], s/2, MPI_INT, 0, MPI_COMM_WORLD);
	
	//pmerge(arrA, arrB, s/2, s/2,myarray); 

	mergeSort(myarray, 0, s);

	//shutdown procedures

	delete[] myarray;
    delete[] arrA;
	delete[] arrB;
	MPI_Finalize();
	return 0;   
}
int Rank(int * a, int first, int last, int valToFind){
	
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
    /*
    for(int i =0; i<indexMerge; i++){
        a[i] = output[i];
    }*/

}

void mergeSort(int * a,int first,int last){
    int m = (first+last)/2;
    int * b = &a[m];
	
	if(m<4) {
        if(first>= m)
            return;
        mergeSort(a, first, m);
        mergeSort(b, first, last-m);
		int* outputarr = new int[last]; 
		smerge(a, b, m, last-m,outputarr);
		for(int i = 0; i< last; i++){
			a[i] = outputarr[i];
		}
		return;
	}
	
    //cout << "let it end " <<first <<endl;
	if( first <last){
		mergeSort(a,first,m);
    	mergeSort(b,first,last-m);
    	int* outputarr = new int[last]; 

    	pmerge(a, b, m, last-m,outputarr);
    //smerge(a, b, m, last-m,outputarr);
		for(int i = 0; i< last; i++){
        	a[i] = outputarr[i];
    	}
	}

}

void pmerge(int *a, int * b, int lasta, int lastb, int * output){
	int currProc = my_rank;
    int numProc = p;
    MPI_Status status;
    

    //srank size 
    int sizeSRank =ceil(((double)lasta)/((double)(log2(lasta))));;
    cout << "Size SRank: " << sizeSRank <<endl;
    //Calculate SRank
    int * A_Select = new int[sizeSRank];
    int * B_Select = new int[sizeSRank];
    int *localSRankA = new int[sizeSRank]; 
    int *localSRankB = new int[sizeSRank]; 
    int * SRankA = new int[sizeSRank];
    int * SRankB = new int[sizeSRank];
    int * localA_Select = new int[sizeSRank];
    int * localB_Select = new int[sizeSRank];
    int * rankA = new int[sizeSRank];
    int * rankB = new int[sizeSRank];

    cout << "Current Processor: " << my_rank << "\na list: ";
    for(int i =0; i <lasta; i++){
        cout << a[i] << " ";
    }
    cout << "\nb list: " ;
    for(int i =0; i <lastb; i++){
        cout << b[i] << " ";
    }
    cout << endl;
     int lm = log2(lasta);
     int lmo = log2(lasta);
    for(int i = 0; i< sizeSRank;i++){
        A_Select[i] = 0;
        B_Select[i] = 0;
        localSRankA[i] =0;
        localSRankB[i] = 0; 
        SRankA[i] = 0;
        SRankB[i] = 0;
        localA_Select[i] = 0;
        localB_Select[i] = 0;
        rankA[i] =0;
        rankB[i] = 0;
    }
    for(int i = 0; i< sizeSRank;i++){
        rankA[i] += lmo;
        rankB[i] += lmo;
        lmo +=log2(lasta);
    }
   
    
    if(my_rank< sizeSRank){
        if(my_rank == 0){
            localA_Select[my_rank] = a[0];
            localB_Select[my_rank] = b[0];
        }
        else{
            localA_Select[my_rank] = a[my_rank*lm];
            localB_Select[my_rank] = b[my_rank*lm];
        }
        localSRankA[my_rank] = Rank(b, 0,lastb, localA_Select[my_rank]);
        localSRankB[my_rank] = Rank(a, 0, lasta, localB_Select[my_rank]);
        //rankA[my_rank] = Rank(a, 0, lasta, localA_Select[my_rank]);
        //rankB[my_rank] = Rank(b, 0, lastb, localB_Select[my_rank]);
    }
    MPI_Allgather(&localA_Select[my_rank], 1, MPI_INT, &A_Select[0], 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&localB_Select[my_rank], 1, MPI_INT, &B_Select[0], 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&localSRankA[my_rank], 1, MPI_INT, &SRankA[0], 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&localSRankB[my_rank], 1, MPI_INT, &SRankB[0], 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&rankA[my_rank], 1, MPI_INT, &rankA[0], 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&rankB[my_rank], 1, MPI_INT, &rankB[0], 1, MPI_INT, MPI_COMM_WORLD);

    if(my_rank==0){
    cout << "Current Processor: " << my_rank <<endl;
    cout <<"ASelect: ";
    for(int i = 0;i<sizeSRank;i++){
        printf("%d ", A_Select[i]);
    }
    cout <<"\nbSelect: ";
    for(int i = 0;i<sizeSRank;i++){
        printf("%d ", B_Select[i]);
    }
    cout <<"\nSRankA: ";
    for(int i = 0;i<sizeSRank;i++){
        printf("%d ", SRankA[i]);
    }
    cout <<"\nSRankB: ";
    for(int i = 0;i<sizeSRank;i++){
        printf("%d ", SRankB[i]);
    }
    cout <<endl;
     cout <<"\nRankB: ";
    for(int i = 0;i<sizeSRank;i++){
        printf("%d ", rankB[i]);
    }
    cout <<endl;
     cout <<"\nRanka: ";
    for(int i = 0;i<sizeSRank;i++){
        printf("%d ", rankA[i]);
    }
    cout <<endl;
    }
/*
    int *WIN = new int[sizeSRank*2];
    int l, m, n =0;
    cout << "WIN: ";
    while(n!= sizeSRank*2){
        if(SRankA[l]<= SRankB[m]){
            WIN[n++] = SRankA[l++];
        }
        else{
            WIN[n++] = SRankB[m++];
        }
        cout << WIN[n-1] << " ";
    }
    cout << endl;*/
    //shapes part
    //int *WIN = new int[lasta+lastb];
    /*for(int i = 0; i < lasta+lastb; i++){
        WIN[i] = 0;
    }*/

    //if Win has to be overall size 
/*
    int *WIN = new int[s];
    for(int i = 0; i < s; i++){
        WIN[i] = 0;
    }
    int selectvalA=0;
	int selectvalB=0;
	int j=0;
	for(int i=0; i<sizeSRank; i++){
		selectvalA = SRankA[i];
		selectvalB = SRankB[i];
		WIN[selectvalA+j] = a[j];
		WIN[selectvalB+j] = b[j];
		j += ceil((double)(log2(lasta)));
	}
    cout << "WIN: ";
    for(int i =0; i <s; i++){
        cout << WIN[i] << " ";
    }
    cout << endl;
*/

    //if win is just for a +b 
    int *WIN = new int[lasta+lastb];
    for(int i = 0; i < lasta+lastb; i++){
        WIN[i] = 0;
    }
/*
    int selectvalA=SRankA[0];
	int selectvalB=SRankB[0];
    int storedValA =a[0];
    int storedValB =b[0];
	int j=0;
    int l,m =1;
    int temp = 0;
    int tempIndex =0;
    while(j< sizeSRank){
        if(temp == 0){
            if(WIN[selectvalA] == 0){
                WIN[selectvalA] = storedValA;
                selectvalA = SRankA[l];
                storedValA = a[l++];
                j++;
            }
            else{
                if(WIN[selectvalA] > storedValA){
                    temp = WIN[selectvalA];
                    tempIndex = selectvalA;
                    WIN[selectvalA] = storedValA;
                    selectvalA = SRankA[l];
                    storedValA = a[l++];
                }
                else{
                    temp = storedValA;
                    tempIndex = selectvalA;
                    selectvalA = SRankA[l];
                    storedValA = a[l++];
                }
            }
            if(WIN[selectvalB] == 0){
                WIN[selectvalB] = storedValB;
                selectvalB = SRankB[m];
                storedValB = b[m++];
                j++;
            }
            else{
                if(WIN[selectvalB] > storedValB){
                    temp = WIN[selectvalB];
                    tempIndex = selectvalB;
                    WIN[selectvalB] = storedValB;
                    selectvalB = SRankB[m];
                    storedValB = b[m++];
                }
                else{
                    temp = storedValB;
                    tempIndex = selectvalB;
                    selectvalB = SRankB[m];
                    storedValB = b[m++];
                }
            }
        }
        else{
            if(WIN[tempIndex+1] == 0){
                WIN[tempIndex+1]  = temp;
                temp = 0;
                tempIndex =0;
                j++;
            }
            else{
                if(WIN[tempIndex+1] < temp){
                    tempIndex++;
                }
                else{
                    int temp2 = WIN[tempIndex+1];
                    WIN[tempIndex+1] = temp;
                    temp =temp2;
                    tempIndex++;
                }
            }
        }
        

    }
*/

	/*while(j < sizeSRank){
        if(WIN[selectvalA] == 0){
            WIN[selectvalA] = a[l];
            selectvalA = SRankA[l];
            l++;
            j++;
        }
        else{
            if(WIN[selectvalA] > )
        }
	}
    for(int i=1; i<sizeSRank; i++){
        if(WIN[selectvalB] ==0){
		    WIN[selectvalB] = b[m];
		    selectvalB = SRankB[m];
            m++;
        }
        else{

        }
		
		j += ceil((double)(log2(lasta)));
	}*/
    cout << "WIN: ";
    for(int i =0; i <lasta+lastb; i++){
        cout << WIN[i] << " ";
    }
    cout << endl;


    //shapes trial 2
    int shapesPerP = ceil(double(sizeSRank*2)/double(p));
    cout << "Shapes Per P: " << shapesPerP << endl;
    int * alowerBounds = new int [sizeSRank*shapesPerP];
    int * aupperBounds = new int [sizeSRank*shapesPerP];
    int * blowerBounds = new int [sizeSRank*shapesPerP];
    int * bupperBounds = new int [sizeSRank*shapesPerP];
    int auCounter =0;
    int alCounter =0;
    int * sizeB = new int[sizeSRank*shapesPerP];
    int * sizeA = new int[sizeSRank*shapesPerP];
    int w = 0;
    int endProc =my_rank*lm;
    int indexShape =0;
    int srankbcounter =0;
    int arankcounter =0;
    int srankacounter =0;
    int brankcounter =0;
	for(int i=0; i<sizeSRank*shapesPerP; i++){
		alowerBounds[i] =0;
	}
	for(int i=0; i<sizeSRank*shapesPerP; i++){
		aupperBounds[i] =0;
	}
    //this sdoesnt work 
    
    {
		while(srankbcounter < sizeSRank && arankcounter < sizeSRank) {//all of a's lower bounds

			if(SRankB[srankbcounter] <rankA[arankcounter]) {
			alowerBounds[alCounter] = SRankB[srankbcounter];
			/*if(my_rank==0){
			cout<< "albounds in first if " << alowerBounds[alCounter] <<endl;
			cout <<"srankb counter " << srankbcounter <<endl;
			


			
			cout << my_rank << " rank is srankb lower bound: "<<alCounter<<endl;
			}*/
			srankbcounter++;
			alCounter++;
			}//srankB's actual value compared to ranka's actual value
			else if(SRankB[srankbcounter] == rankA[arankcounter])  {
				if(B_Select[srankbcounter] < A_Select[arankcounter]){

                alowerBounds[alCounter]=rankA[arankcounter];
				SrankB[srankbcounter]++;
				/*if(my_rank==0){
					cout<< "albounds in second  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;
				}*/
	            arankcounter++;
				alCounter++;
                }
            
				
                
				if(B_Select[srankbcounter] > A_Select[arankcounter]){

                alowerBounds[alCounter]=SRankB[srankbcounter];
				rankA[arankcounter]++;
				/*if(my_rank==0){
					cout<< "albounds in second  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;
				}*/

                
            	srankbcounter++;
				alCounter++;
                    }
				}
                
                else if(SRankB[srankbcounter] >rankA[arankcounter])  {
				alowerBounds[alCounter]=rankA[arankcounter];
				
				/*if(my_rank==0){
					cout<< "albounds in second  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;


				}*/
				arankcounter++;
				alCounter++;
				}
			}
			{
			while(arankcounter < sizeSRank){
				alowerBounds[alCounter]=rankA[arankcounter];
				

				/*if(my_rank==0){
					cout<< "albounds in third  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;


				}*/
				arankcounter++;
				alCounter++;
				}

			}
			{
			while(srankbcounter < sizeSRank){
				alowerBounds[alCounter]=SRankB[srankbcounter];
			/*if(my_rank==0){
				cout<< "albounds in fourth  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is srankb lower bound: "<<alCounter<<endl;
				cout <<"srankb counter " << srankbcounter <<endl;

				}*/


				srankbcounter++;
				alCounter++;
				}

			}
			if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
				
				if(alowerBounds[sizeSRank*shapesPerP-1] >= lasta){

				alowerBounds[sizeSRank*shapesPerP-1] = lasta;
			}
				
				cout  <<" A lowerbounds: ";
			for(int i=0; i<sizeSRank*shapesPerP; i++){
				cout<< alowerBounds[i] << " ";
			}
			
			cout<<endl;
			MPI_Bcast(&alowerBounds[0], sizeSRank*shapesPerP, MPI_INT, 0, MPI_COMM_WORLD);
			}

	}
//b lower bound ------------------------------------------------
	{
		alCounter=0;
		while(srankacounter < sizeSRank && brankcounter < sizeSRank) {//all of a's lower bounds
			if(SRankA[srankacounter] <rankB[brankcounter]) {
				blowerBounds[alCounter] = SRankA[srankacounter];
				/*if(my_rank==0){
				cout<< "albounds in first if " << alowerBounds[alCounter] <<endl;
				cout <<"srankb counter " << srankbcounter <<endl;
				cout << my_rank << " rank is srankb lower bound: "<<alCounter<<endl;
				}*/
				srankacounter++;
				alCounter++;

			}
           //srankB's actual value compared to ranka's actual value
			else if(rankB[brankcounter] == SRankA[srankacounter])  {
				if(B_Select[brankcounter] < A_Select[srankacounter]){

                alowerBounds[alCounter]=SRankA[srankacounter];
				rankB[brankcounter]++;
				/*if(my_rank==0){
					cout<< "albounds in second  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;
				}*/
	            brankcounter++;
				alCounter++;
                }
            
				
                
				if(B_Select[srankbcounter] > A_Select[arankcounter]){

                alowerBounds[alCounter]=SRankB[srankbcounter];
				rankA[arankcounter]++;
				/*if(my_rank==0){
					cout<< "albounds in second  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;
				}*/

                
            	srankbcounter++;
				alCounter++;
                    }
				}
			else {
				blowerBounds[alCounter]=rankB[brankcounter];
				/*if(my_rank==0){
					cout<< "albounds in second  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;
				}*/
				brankcounter++;
				alCounter++;
			}
		}
		
		while(brankcounter < sizeSRank){
			blowerBounds[alCounter]=rankB[brankcounter];
				/*if(my_rank==0){
					cout<< "albounds in third  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is arank lower bound: "<<alCounter<<endl;
				cout <<"arank counter " << arankcounter <<endl;


				}*/
				brankcounter++;
				alCounter++;
		}	
		while(srankacounter < sizeSRank){
			blowerBounds[alCounter]=SRankA[srankacounter];
			/*if(my_rank==0){
				cout<< "albounds in fourth  " << alowerBounds[alCounter] <<endl;
				cout << my_rank << " rank is srankb lower bound: "<<alCounter<<endl;
				cout <<"srankb counter " << srankbcounter <<endl;

			}*/
			srankacounter++;
			alCounter++;
		}

			
		if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
			
				if(blowerBounds[sizeSRank*shapesPerP-1] >= lastb){

				blowerBounds[sizeSRank*shapesPerP-1] = lastb;
			}
			
			cout  <<" B lowerbounds: ";
			for(int i=0; i<sizeSRank*shapesPerP; i++){
				cout<< blowerBounds[i] << " ";
			}
			cout<<endl;
			MPI_Bcast(&blowerBounds[0], sizeSRank*shapesPerP, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	//Upper bounds -- A 
	{
		int indexAlBounds = 0;
		srankbcounter= 0;
		arankcounter = 0;
		while(indexAlBounds != sizeSRank*shapesPerP){
			if(alowerBounds[indexAlBounds] ==0){ // do this for if the lower bounds ==0 then upper bound = 0 -- i might be wrong for this
				aupperBounds[auCounter] = 0;
				auCounter++;
				indexAlBounds++;
				//have to increment corrisponding arankbcounter and arank counter 
				if(rankA[arankcounter] == 0){
					arankcounter++;
				}
				else{
					srankbcounter++;
				}
			}
			else{ //do this if alBounds[i] != 0
				if(alowerBounds[indexAlBounds] == SRankB[srankbcounter]){
					aupperBounds[auCounter] = alowerBounds[indexAlBounds-1]+1;
					auCounter++;
					indexAlBounds++;
					srankbcounter++;
				}
				else{// if(alowerBounds[indexAlBounds] == rankA[arankcounter]){
					aupperBounds[auCounter] = alowerBounds[indexAlBounds-1]+1;
					auCounter++;
					indexAlBounds++;
					arankcounter++;
				}
				
			}
			
		}
		if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
				if(aupperBounds[sizeSRank*shapesPerP-1] >= lasta){

				aupperBounds[sizeSRank*shapesPerP-1] = lasta;
			}
			cout  <<" A upperbounds: ";
			for(int i=0; i<sizeSRank*shapesPerP; i++){
				cout<< aupperBounds[i] << " ";
			}
			cout<<endl;
			MPI_Bcast(&aupperBounds[0], sizeSRank*shapesPerP, MPI_INT, 0, MPI_COMM_WORLD);
		}
		

	}
	//Upper bounds - B
	{
		int indexBlBounds = 0;
		auCounter =0;
		srankacounter = 0;
		brankcounter = 0;
		while(indexBlBounds != sizeSRank*shapesPerP){
			if(blowerBounds[indexBlBounds] ==0 ){ // do this for if the lower bounds ==0 then upper bound = 0 -- i might be wrong for this
				//have to increment corrisponding arankbcounter and arank counter 
				if(rankB[brankcounter] == 0){
					bupperBounds[auCounter] = 0;
					brankcounter++;
				}
				else if(SRankA[srankacounter] ==0){
					bupperBounds[auCounter] = 0;
					srankacounter++;
				}
				auCounter++;
				indexBlBounds++;
			}
			else{ //do this if alBounds[i] != 0
				if(blowerBounds[indexBlBounds] == SRankA[srankacounter]){
					bupperBounds[auCounter] = blowerBounds[indexBlBounds-1]+1;
					auCounter++;
					indexBlBounds++;
					srankacounter++;
				}
				else{// if(alowerBounds[indexAlBounds] == rankA[arankcounter]){
					bupperBounds[auCounter] = blowerBounds[indexBlBounds-1]+1;
					auCounter++;
					indexBlBounds++;
					brankcounter++;
				}
				
			}
			
		}
		if(my_rank==0){//need to look at double zeros and compare them to see which is bigger for the shape
				if(bupperBounds[sizeSRank*shapesPerP-1] >= lastb){

				bupperBounds[sizeSRank*shapesPerP-1] = lastb;
			}
			cout  <<" B upperbounds: ";
			for(int i=0; i<sizeSRank*shapesPerP; i++){
				cout<< bupperBounds[i] << " ";
			}
			cout<<endl;
	


			MPI_Bcast(&bupperBounds[0], sizeSRank*shapesPerP, MPI_INT, 0, MPI_COMM_WORLD);
		}
		
		

	}
    
    /*
    for(int i =0; i<shapesPerP; i++){
        if(i == shapesPerP -1){
            alowerBounds[alCounter] = ((my_rank+shapesPerP)*lm)-1;
            if(SRankB[my_rank+shapesPerP] <= alowerBounds[alCounter]){
                aupperBounds[auCounter] = SRankB[my_rank+shapesPerP];
                auCounter++;
            }
            else{
                aupperBounds[auCounter] = my_rank+i;
            }
        }
    }*/
   // MPI_Allgather(&alowerBounds[my_rank], shapesPerP, MPI_INT, &alowerBounds[0], shapesPerP, MPI_INT, MPI_COMM_WORLD);
    //MPI_Allgather( const void* sendbuf , MPI_Count sendcount , MPI_Datatype sendtype , void* recvbuf , MPI_Count recvcount , MPI_Datatype recvtype , MPI_Comm comm);
    //MPI_Allreduce(&alowerBounds[my_rank], &alowerBounds[my_rank], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
    

}