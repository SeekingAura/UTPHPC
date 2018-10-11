#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <chrono>
#include <omp.h>
using namespace std;



class matrix{
	public:
	int *M1, *M2, *MResult; /* matrices (M1,M2) */
	size_t M1row=0,M1col=0, M2row=0, M2col=0;
	//constructor
	matrix(){

	}

	matrix(const char *fileName1, bool useMR){
		FILE *f1=NULL; /* file pointers */
		//this->buildMatrix();
		//printf("constructor\n");
		f1=this->openFile(fileName1);

		this->M1=this->buildMatrix(f1, this->M1row, this->M1col);
		this->getData(f1, this->M1, this->M1row*this->M1col);

		this->M2=this->buildMatrix(f1, this->M2row, this->M2col);
		this->getData(f1, this->M2, this->M2row*this->M2col);
		this->MResult=new int[this->M1row*this->M2col];
		
		
		
		fclose(f1);
		//f2=this->openFile(fileName2);


		// fscanf(f,"%i",this->M2row); /* %zu zx is size_t */
		// fscanf(f,"%i",this->M2col);	

	}
	//destructor
	~matrix(){
		delete [] this->M1;
		delete [] this->M2;
		delete [] this->MResult;
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

	int * buildMatrix(FILE *f, size_t &rows, size_t &columns){
		/* build a matrix M (get memory) */
		fscanf(f,"%i",&rows); /* %zu zx is size_t */
		fscanf(f,"%i",&columns);
		//fgetc(f);  /* skipping nasty character, on this case new line */
		int *M;
		M=new int[rows*columns];
		if(rows <= 0 || columns <= 0){
			printf("should be size positive and upper than cero \n");
			exit(1);
		}
		return M;
	}

	void getData(FILE *f, int *M, size_t len){
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
		char data[10]="", ch = ' ';
		size_t posData = 0, Mindex = 0;
		while(len>Mindex){
			ch = fgetc(f); /*get char and char in file f */
			if(Mindex==0 && ch == '\n'){//skip nasty chracter
				continue;
			}
			if(ch == ',' || ch == '\n'){
				data[posData] = '\0'; /* char end */
				M[Mindex] = atoi(data); /*convert string to int */
				posData = 0;
				strcpy(data, ""); /* take memory for the next data */
				Mindex++;
			}else{
				data[posData] = ch;
				posData++;
			}

		}
	}

	void buildMatrixTemp(){//reserve memory for Workers data
		this->M1= new int[this->M1row*this->M1col];
		this->M2= new int[this->M2row*this->M2col];
		this->MResult= new int[this->M1row*this->M2col];

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
					fprintf(f,"%i\n",this->MResult[i*this->M2col + j]);
				}
				else {
					fprintf(f,"%i,",this->MResult[i*this->M2col + j]);
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

	void mulSecuential(int *MResult){
		/*
			This function will multiply two matrices (M1,M2)
			M1 -> Matrix1, M2 -> Matrix2, M1row -> Matrix1 rows, M1col -> Matrix1
			columns, M2row -> Matrix2 rows, M2col -> Matrix2 columns
		*/

		for(size_t i=0; i<this->M1row; i++){//every row m1
			for(size_t j=0; j<this->M2col; j++){//every column m2
				int data = 0;
				for(size_t k=0; k<this->M1col; k++){//take in row per value column
					data = this->M1[i*this->M1col+k] * this->M2[k*this->M2col+j] + data;
				}
				MResult[i*this->M1col+j] = data;
			}
		}
	}

	void mulParallelRow(int *MResult, int *M1, int *M2, int M1row, int M2col){
		/*
			This function will multiply two matrices (M1,M2)
			M1 -> Matrix1, M2 -> Matrix2, M1row -> Matrix1 rows, M1col -> Matrix1
			columns, M2row -> Matrix2 rows, M2col -> Matrix2 columns
		*/
		size_t i=0,j=0,k=0, chunk=this->M1row/8;
		int numThreads = 0;
		if(chunk==0){
			chunk=1;
		}
		#pragma omp parallel private(i,j,k,numThreads) shared(MResult,chunk) //num_threads(8)
		{
			#pragma omp for schedule(dynamic,chunk)
				for(i=0; i<M1row; i++){//every row m1
					for(j=0; j<M2col; j++){//every column m2
						int data = 0;
						for(k=0; k<M2col; k++){//take in row per value column
							data = M1[i*M2col+k] * M2[k*M2col+j] + data;
						}
						MResult[i*M2col+j] = data;
					}
				}
		}
		
	}



	void printMatrix(int *M, size_t row, size_t col){
		/* getting pos of matrix
			->> rowi
			0 | 1 | 2
			3 | 4 | 5
			columni
			rows -> 2
			columns -> 3
			rowi*columns+columni=pos
		*/
		for(size_t i=0; i<row; i++){
			for(size_t j=0; j<col; j++){
				if(j==col-1){
					printf("%i", M[i*col+j]);
				}else{
					printf("%i | ", M[i*col+j]);
				}	
			}
			printf("\n");
		}
	}

	void printOperators(){
		printf("matrix M1 -> \n");
		this->printMatrix(this->M1, this->M1row, this->M1col);

		printf("matrix M2 -> \n");
		this->printMatrix(this->M2, this->M2row, this->M2col);
	}

	void printResult(){
		printf("matrix MResult -> \n");
		this->printMatrix(this->MResult, this->M1row, this->M2col);
	}

	void printResult(int *M){
		printf("matrix MResult -> \n");
		this->printMatrix(M, this->M1row, this->M2col);
	}

};




void writeTime(double elapsed, size_t len){
		/*
			Wite the time on timesc++Parallel.txt file
			M -> Matrix, Mrow -> Matrix rows, Mcol -> Matrix columns
		*/
		FILE *f = fopen("timesc++Parallel.txt","a+");//write at end of file and set result, append
		fprintf(f,"%i	%.9f\n", len, elapsed);
		fclose(f);
}

void executeHeader(){
	matrix opMatrix(argv[1], 0);
	int stepPart=opMatrix.M1row/(p-1);
	
	for(int nodeWorkerId=1, startPart=0, endPart=opMatrix.M1row/(p-1);nodeWorkerId<=p;nodeWorkerId++, endPart+=stepPart){
		
		if(nodeWorkerId==p){
			endPart=opMatrix.M1row;
		}
		MPI_Send(endPart-startPart, 1, MPI_INT, nodeWorkerId, MSGTAG, MPI_COMM_WORLD);
		MPI_Send(opMatrix.M1col, 1, MPI_INT, nodeWorkerId, MSGTAG, MPI_COMM_WORLD);
		Mtemp= new int[(endPart-startPart)*opMatrix.M1col]
		for(int numRow=startPart numPos=0; numRow<endPart; numRow++){//Fill part of MatrixTemp (operator rows)
			for(int numCol=0; numCol<opMatrix.M1col; numCol++, numPos++){
				Mtemp[numPos]=opMatrix.M1[numRow*opMatrix.M1row+numCol]
			}
		}
		MPI_Send(opMatrix.M2row, 1, MPI_INT, nodeWorkerId, MSGTAG, MPI_COMM_WORLD);
		MPI_Send(opMatrix.M2col, 1, MPI_INT, nodeWorkerId, MSGTAG, MPI_COMM_WORLD);
		MPI_Send(MTemp, endPart-startPart,MPI_INT, nodeWorkerId, MSGTAG, MPI_COMM_WORLD);
		MPI_Send(opMatrix.M2, opMatrix.M2row*opMatrix.M2col, MPI_INT, nodeWorkerId, MSGTAG, MPI_COMM_WORLD);
		startPart=endPart;
		delete [] Mtemp;
	}
	
	for(int nodeWorkerId=1, startPart=0, endPart=opMatrix.M1row/(p-1);nodeWorkerId<=p;nodeWorkerId++, endPart+=stepPart){
		if(nodeWorkerId==p){
			endPart=opMatrix.M1row;
		}
		MPI_Recv(MTemp, endPart-startPart, MPI_INT, nodeWorkerId, MSGTAG, MPI_COMM_WORLD, &status);
		for(int numRow=startPart numPos=0; numRow<endPart; numRow++){//Fill part of MatrixTemp (operator rows) into MResult
			for(int numCol=0; numCol<opMatrix.M1col; numCol++, numPos++){
				opMatrix.MResult[numRow*opMatrix.M1row+numCol]=Mtemp[numPos];
			}
		}
	}

	opMatrix.printResult();
}

void executeWorker(){
	int NodeHeaderId=0;
	matrix opMatrix();
	// M1 info
	MPI_Recv(opMatrix.M1row, 1, MPI_INT, NodeHeaderId, MSGTAG, MPI_COMM_WORLD, &status);
	MPI_Recv(opMatrix.M1col, 1, MPI_INT, NodeHeaderId, MSGTAG, MPI_COMM_WORLD, &status);

	//M2 info
	MPI_Recv(opMatrix.M2row, 1, MPI_INT, NodeHeaderId, MSGTAG, MPI_COMM_WORLD, &status);
	MPI_Recv(opMatrix.M2col, 1, MPI_INT, NodeHeaderId, MSGTAG, MPI_COMM_WORLD, &status);

	opMatrix.buildMatrixTemp();
	MPI_Recv(opMatrix.M1, opMatrix.M1row*opMatrix.M1col, MPI_INT, NodeHeaderId, MSGTAG, MPI_COMM_WORLD, &status);
	MPI_Recv(opMatrix.M2, opMatrix.M2row*opMatrix.M2col, MPI_INT, NodeHeaderId, MSGTAG, MPI_COMM_WORLD, &status);

	opMatrix.mulParallelRow(opMatrix.MResult, opMatrix.M1, opMatrix.M2, opMatrix.M1row, opMatrix.M2col);

	MPI_Send(opMatrix.MResult, opMatrix.M1row*opMatrix.M2col, MPI_INT, NodeHeaderId, MSGTAG, MPI_COMM_WORLD);
}

void writeTime(float elapsed, int valueK, int len){
	FILE *f = fopen("timesc++MPI.txt","a+");//write at end of file and set result, append
	fprintf(f,"%ld	%.9f\n", valueK, len, elapsed);
	fclose(f);
}

int main(int argc, char const *argv[]) {
	if(argc != 2){
		printf("There should be 2 arguments!\n");
		exit(1);
	}
	int p_id;
  	int p;
	int MSGTAG=0;

	MPI_Status status;
	MPI_Init ( &argc, &argv );
	MPI_Comm_rank ( MPI_COMM_WORLD, &p_id );//identifica el número de equipo que está corriendo
	MPI_Comm_size ( MPI_COMM_WORLD, &p );//identifica el total de equipos que se usarán
	auto startTime=std::chrono::high_resolution_clock::now();
	if(p_id==0){//Header Part
		executeHeader();
	}else{//Workers part
		executeWorker()
	}
	auto endTime=std::chrono::high_resolution_clock::now();
	writeTime(elapsed.count(), n, nsteps);

	double startTime = omp_get_wtime();
	opMatrix.mulParallelRow(opMatrix.MResult, opMatrix.M1, opMatrix.M2, opMatrix.M1row, opMatrix.M2col);
	double endTime = omp_get_wtime();
	double elapsed = endTime - startTime;
	writeTime(elapsed, opMatrix.M1row);
	
	// opMatrix.printOperators();
	// opMatrix.printResult();
	return 0;
}