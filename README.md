This project involves implementing a parallel version of the mergesort algorithm in C++ using MPI (Message Passing Interface). The goal is to create a "clone" of the traditional mergesort, leveraging parallel processing to improve performance. The implementation includes writing and testing several functions that divide, sort, and merge arrays, while simulating shared memory between processes.

Key Functions:
*mergesort(int a, int first, int last)**: The primary sorting function which will divide and conquer the array into smaller sub-arrays.
smerge(int a, int b, int lasta, int lastb, int* output = NULL)**: A merging function that combines two sorted sub-arrays into a single sorted array.
*rank(int a, int first, int last, int valToFind)**: Determines the rank (position) of a value in a sorted array, used to facilitate the merging process.
pmerge(int a, int b, int lasta, int lastb, int* output = NULL)**: The parallel merge function, which merges two sorted arrays using multiple processors, simulating parallelism by striping work across processes.

Implementation Requirements:
Dynamically allocate an array based on user input, and populate it with random numbers (e.g., between 0 and 500).
Use MPI to simulate shared memory, where process 0 creates the initial array and broadcasts it to other processes.
Implement and test the smerge and rank functions independently before moving on to the parallel merge.
The parallel mergesort must handle arrays up to size 128, with at least two test runs demonstrating correct sorting for an array of that size.

Phases of Parallel Merging (pmerge):
Calculate SRANKA and SRANKB: Distribute the work of computing sample ranks across processors.
Share Ranks: Ensure that all processes have the correct rank values.
Place Sampled Elements: Insert sampled rank elements into their correct positions in an output array.
Shape Game: Divide the work based on shapes derived from the rank information, processed in parallel using smerge.
Final Merge: Use MPI to reduce and gather the results, ensuring all processes share the same final sorted array.

Testing & Debugging:
Test each function thoroughly, especially when handling parallel tasks and recursion in the merge.
Debug the solution by working with smaller arrays and using fixed seeds for reproducibility.
