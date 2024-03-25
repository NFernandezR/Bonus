#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
// Implement heapify function to implement heap sort.
// Heapify subtree at the i-th index value, for a heap of size n.
void heapify(int arr[], int n, int i)
{
	// Initialize the biggest value to be that of the root index of the subtree.
	int biggest = i;

	// Initialize the variable storing the left child node's value.
	int left = (2 * i) + 1;
	// Initialize the variable storing the right child node's value.
	int right = (2 * i) + 2;

	// If the left node's value is greater than the current biggest...
	if (left < n && arr[left] > arr[biggest])
	{
		// ... update the biggest value to be the left node's value.
		biggest = left;
	}

	// If the right node's value is greater than the current biggest...
	if (right < n && arr[right] > arr[biggest])
	{
		// ... update the biggest value to be the right node's value.
		biggest = right;
	}

	// If the root value is not the biggest value, swap the values, then continue heapifying.
	if (biggest != i)
	{
		// Swap the root (i-th) and biggest values.
		int temp = arr[i];
		arr[i] = arr[biggest];
		arr[biggest] = temp;

		// Recursively heapify subtree
		heapify(arr, n, biggest);
	}
}


void heapSort(int arr[], int n)
{
	// Heapify the array to obtain the maximum heap.
	for (int i = n / 2 - 1; i >= 0; i--)
	{
		heapify(arr, n, i);
	}

	// Sort the heap
	for (int j = n - 1; j >= 0; j--)
	{
		// Swap the root and first elements.
		int temp = arr[0];
		arr[0] = arr[j];
		arr[j] = temp;

		// Heapify the root (j-th) item to find the biggest value at the root.
		heapify(arr, j, 0);
	}
}

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r)
{
	// Base Case: Return if the array or subarray is sized for only a single element.
	// Recursive Case: Subdivide the array into halves until the base case is reached, then sort those subarrays,
	// backtracking until the entire array is sorted.
	if (l < r)
    {
		// Compute the middlemost index position of the array.
    	int m = (l + r) / 2;

        // Recursively call mergeSort to sort first the left, then the right halves of the array.
        mergeSort(pData, l, m);
        mergeSort(pData, m + 1, r);
        
        // Merge the sorted left and right halves into a sorted, whole array.
    	// Array Index Variables
    	int i, j, k;
    	// Compute the size of temporary subarray, L. Store as an int, left.
    	int left = m - l + 1;
    	// Size of temporary subarray, R. Store as an int, right.
    	int right = r - m;

    	// Allocate memory for the temporary subarrays, L and R.
    	int *L = Alloc(sizeof(int)*left);
    	int *R = Alloc(sizeof(int)*right);

    	// Copy the corresponding data from the original array over to the subarrays.
    	for (i = 0; i < left; i++)
    	{
			// Fill each i-th entry of the left-half subarray, L, with the data at the "l + i"-th position of the original array.
        	L[i] = pData[l + i];
    	}
    	for (j = 0; j < right; j++)
    	{
			// Fill each j-th entry of the right-half subarray, R, with the data at the "m + 1 + j"-th position of the original array.
        	R[j] = pData[m + 1 + j];
    	}

    	// Fill the original array with the entries from the temp arrays in sorted order
    	// Initialize array index variables
    	i = 0; j = 0; k = l;

    	// While BOTH the left- and right-half subarray have not been COMPLETELY iterated over...
		while (i < left && j < right)
    	{
			// If the i-th entry of L is less than or equal to the j-th entry of R...
        	if (L[i] <= R[j])
        	{
				// The k-th entry of the original array takes on the data stored in the i-th entry of L.
            	pData[k] = L[i];
				// Move along to the next element in the left subarray, loop ends if no next element exists.
            	i++;
        	}

			// Else, the i-th entry of L is greater than the j-th entry of R, so...
        	else
        	{
				// The k-th entry of the original array takes on the data stored in the j-th entry of R.
            	pData[k] = R[j];
				// Move along to the next element in the right subarray, loop ends if no next element exists.
            	j++;
        	}
			// Move along to the next element of the original array.
        	k++;
    	}

    	// It is possible that one subarray was larger in size than the other, so to address these leftover entries...
		// While the left subarray has not been completely iterated over...
    	while (i < left)
    	{
			// Assign the i-th entry of L to the k-th entry of the original array.
        	pData[k] = L[i];
			// Move along to the next elements of the left subarray and original array.
        	i++;
        	k++;
    	}
		// While the right subarray has not been completely iterated over...
    	while (j < right) {
			// Assign the j-th entry of R to the k-th entry of the original array.
        	pData[k] = R[j];
			// Move along to the next elements of the right subarray and original array.
        	j++;
        	k++;
    	}

    	// Deallocate the memory reserved for the temporary subarrays now that they are no longer needed.
    	DeAlloc(L);
    	DeAlloc(R);
    }
}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n)
{
	int i, item, j;
	for (i = 1; i < n; i++)
	{
		item = pData[i];
		// move
		for (j = i - 1; j >= 0; j--)
		{
			if (pData[j] > item)
			{
				pData[j + 1] = pData[j];
			}
			else
			{
				break;
			}
		}
		pData[j + 1] = item;
	}
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{
	int i, j, temp;
	for (i = 0; i < (n - 1); i++)
	{
		for (j = 0; j < (n - i - 1); j++)
		{
			if (pData[j] > pData[j + 1])
			{
				temp = pData[j];
				pData[j] = pData[j + 1];
				pData[j + 1] = temp;
			}
		}
	}
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n)
{
	// Declare loop variables i and j, minimum index variable mindex, and a temp variable for swapping values.
	int i, j, mindex, temp;
	for (i = 0; i < (n - 1); i++)
	{
		// At the start of every iteration, initialize the minimum index to equal the outter loop iteration number.
		mindex = i;

		for (j = (i + 1); j < n; j++)
		{
			// If the value at the jth index is less than the value at the minimum index...
			if (pData[j] < pData[mindex])
			{
				// ... Assign the minimum index to now be the jth index.
				mindex = j;
			}
		}

		// Swap the value at the ith index with the value at the minimum index.
		temp = pData[i];
		pData[i] = pData[mindex];
		pData[mindex] = temp;
	}
}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
	FILE* inFile = fopen(inputFileName,"r");
	int dataSz = 0;
	*ppData = NULL;
	
	if (inFile)
	{
		fscanf(inFile,"%d\n",&dataSz);
		*ppData = (int *)Alloc(sizeof(int) * dataSz);
		// Implement parse data block
	}
	
	return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

                printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
