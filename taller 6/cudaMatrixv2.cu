#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#define rows 1000
//#define cols 1000
using namespace std;

// CUDA kernel. Each thread takes care of one element of c
__global__ void matricesMul(double *m1, double *m2, double *m3){
		// Get our global thread ID
		int ti = blockIdx.y*blockDim.y+threadIdx.y;
		int tj = blockIdx.x*blockDim.x+threadIdx.x;
		// Make sure we do not go out of bounds
		if(ti < rows && tj < cols){
			double data= 0.0;
			for(int k=0;k<rows;k++) data += m1[ti*rows+k] * m2[k*cols+tj];
			m3[ti*rows+tj] = data;
		}
}

FILE * openFile(char const *fileName){
	/* try to open a file */
	FILE *f=NULL;
	f = fopen(fileName,"r");
	if(f == NULL){
		printf("File '%s' doesn't exist!\n",fileName);
		exit(1);
	}
	return f;
}

double * buildMatrix(FILE *f, size_t &rows, size_t &columns){
	/* build a matrix M (get memory) */
	fscanf(f,"%i",&rows); /* %zu zx is size_t */
	fscanf(f,"%i",&columns);
	//fgetc(f);  /* skipping nasty character, on this case new line */
	double *M = (double *)malloc(rows*columns*sizeof(double)); /* reserved memory */
	if(rows <= 0 || columns <= 0){
		printf("should be size positive and upper than cero \n");
		exit(1);
	}
	return M;
}

void getData(FILE *f, double *M, size_t len){
	/* Capture data from plain text file to system memory
	Note: the data files need one end line to get last number
	format of data files 
	...
	A
	B
	#.#,#.#
	#.#,#.#

	...
	A -> Size row
	B -> Size column
	*/
	//sizeof(char)==1
	char data[50]="", ch = ' ';
	size_t posData = 0, Mindex = 0;
	while(len>Mindex){
		ch = fgetc(f); /*get char and char in file f */
		if(Mindex==0 && ch == '\n'){//skip nasty chracter
			continue;
		}
		if(ch == ',' || ch == '\n'){
			data[posData] = '\0'; /* char end */
			M[Mindex] = stod(data); /*convert string to double */
			posData = 0;
			strcpy(data, ""); /* take memory for the next data */
			Mindex++;
		}else{
			data[posData] = ch;
			posData++;
		}

	}
}
void writeResult(){
	/*
		Wite the result on output.txt file
		M -> Matrix, Mrow -> Matrix rows, Mcol -> Matrix columns
	*/
	FILE *f = fopen("output.txt","w+");//clean file and set result
	for(size_t i=0;i<this->M1row;i++){
		for(size_t j=0;j<this->M2col;j++){
			if(j+1 == this->M2col) {//last chracter
				fprintf(f,"%f\n",this->MResult[i*this->M2col + j]);
			}
			else {
				fprintf(f,"%f,",this->MResult[i*this->M2col + j]);
			}
		}
	}
	fclose(f);
}

bool checkMul(){
	if(this->M1col != this->M2row){
		printf("ERROR - Matrices cannot be multiply!"); 
		return 0;//FALSE
	}
	return 1;//TRUE
}

int main( int argc, char* argv[] ){
	if(argc != 2){
		printf("There should be 2 arguments!\n");
		exit(1);
	}
	// Host (CPU) input matrices
	double *h_m1;
	size_t rows_m1, cols_m1;
	double *h_m2;
	size_t rows_m2, cols_m2;
	//Host (CPU) output matrix
	double *h_m3;

	// Device (GPU-Nvidia) input matrices
	double *d_m1;
	double *d_m2;
	//Device (GPU-Nvidia) output matrix
	double *d_m3;

	FILE *f1=NULL; /* file pointers */
	f1=openFile(argv[1]);
    // Allocate memory for each matrix on host
	h_m1=buildMatrix(f1, rows_m1, cols_m1);
	
	getData(f1, h_m1, rows_m1*cols_m1);

	h_m2=buildMatrix(f1, rows_m2, cols_m2);
	getData(f1, h_m2, rows_m2*cols_m2);

	h_m3=(double *)malloc(rows_m1*cols_m2*sizeof(double));
	
	fclose(f1);
	// Size of matrices n²
    size_t n = rows_m1*cols_m2;
    size_t bytes = n*sizeof(double);

	// Allocate memory for each matrix on GPU
	cudaMalloc((void **)&d_m1, bytes);
	cudaMalloc((void **)&d_m2, bytes);
	cudaMalloc((void **)&d_m3, bytes);


	// Copy host matrices to device
	cudaMemcpy( d_m1, h_m1, bytes, cudaMemcpyHostToDevice);
	cudaMemcpy( d_m2, h_m2, bytes, cudaMemcpyHostToDevice);

	// Number of threads in each thread matrix block
	double x = sqrt(1024);
	size_t threadsInX= floor(x);
	size_t threadsInY= threadsInX;
	dim3 dimBlock(threadsInX,threadsInY,1);
	// Number of thread blocks in matrix grid
	size_t gridNum = ceil((double)n/1024);  // needed grid numbers to our problem
	size_t gridR = ceil(sqrt(gridNum)); 		// grid rows	
	size_t gridC = gridR;										// grid cols
	dim3 dimGrid(gridR,gridC,1);

	// Execute the kernel
	matricesMul<<<dimGrid,dimBlock>>>(d_m1, d_m2, d_m3);

	// Copy result m3 matrix back to host
	cudaMemcpy(h_m3, d_m3, bytes, cudaMemcpyDeviceToHost);

	// print every item into m3 matrix
	for(int i=0; i<n; i++){
		double val = h_m3[i];
		printf("final result: %f\n", val);
	}

	// Release device memory
	cudaFree(d_m1);
	cudaFree(d_m2);
	cudaFree(d_m3);

	// Release host memory
	free(h_m1);
	free(h_m2);
	free(h_m3);

	return 0;
}
