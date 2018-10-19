#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <chrono>
#include <omp.h>

using namespace std;

FILE * openFile(char *fileName){
	/* try to open a file */
	FILE *f=NULL;
	f = fopen(fileName,"r");
	if(f == NULL){
		printf("File '%s' doesn't exist!\n",fileName);
		exit(1);
	}
	return f;
}

void getData(FILE *f, int *M, int len){
	/* 
		Capture data from plain text file to system memory
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
	int posData = 0, Mindex = 0;
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

int * buildMatrix(FILE *f, int &rows, int &columns){
	/* build a matrix M (get memory) */
	fscanf(f,"%i",&rows); /* %zu zx is int */
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

void constructMatrix(char *fileName1, int &M1row, int &M1col, int &M2row, int &M2col, int *M1, int *M2){
	FILE *f1=NULL; /* file pointers */
	//this->buildMatrix();
	//printf("constructor\n");
	f1=openFile(fileName1);

	M1=buildMatrix(f1, M1row, M1col);
	getData(f1, M1, M1row*M1col);

	M2=buildMatrix(f1, M2row, M2col);
	getData(f1, M2, M2row*M2col);
	
	
	
	fclose(f1);
	//f2=this->openFile(fileName2);


	// fscanf(f,"%i",this->M2row); /* %zu zx is int */
	// fscanf(f,"%i",this->M2col);	

}


	
void writeResult(int M1row, int M2col, int *MResult){
	/*
		Wite the result on output.txt file
		M -> Matrix, Mrow -> Matrix rows, Mcol -> Matrix columns
	*/
	FILE *f = fopen("output.txt","w+");//clean file and set result
	for(int i=0;i<M1row;i++){
		for(int j=0;j<M2col;j++){
			if(j+1 == M2col) {//last chracter
				fprintf(f,"%i\n", MResult[i*M2col + j]);
			}
			else {
				fprintf(f,"%i,",MResult[i*M2col + j]);
			}
		}
	}
	fclose(f);
}



bool checkMul(int M1col, int M2row){
	if(M1col != M2row){
		printf("ERROR - Matrices cannot be multiply!"); 
		return 0;//FALSE
	}
	return 1;//TRUE
}

void mulSecuential(int M1row, int M1col, int M2row, int M2col, int *M1, int *M2, int *MResult){
	/*
		This function will multiply two matrices (M1,M2)
		M1 -> Matrix1, M2 -> Matrix2, M1row -> Matrix1 rows, M1col -> Matrix1
		columns, M2row -> Matrix2 rows, M2col -> Matrix2 columns
	*/

	for(int i=0; i<M1row; i++){//every row m1
		for(int j=0; j<M2col; j++){//every column m2
			int data = 0;
			for(int k=0; k<M1col; k++){//take in row per value column
				data = M1[i*M1col+k]*M2[k*M2col+j] + data;
			}
			MResult[i*M1col+j] = data;
		}
	}
}

void mulParallelRow(int startRow, int endRow, int M2col, int *M1, int *M2, int *MResult){
	/*
		This function will multiply two matrices (M1,M2)
		M1 -> Matrix1, M2 -> Matrix2, M1row -> Matrix1 rows, M1col -> Matrix1
		columns, M2row -> Matrix2 rows, M2col -> Matrix2 columns
	*/
	int i=0,j=0,k=0, chunk= (endRow - startRow) /8;
	int numThreads = 0;
	if(chunk==0){
		chunk=1;
	}
	#pragma omp parallel private(i,j,k,numThreads) shared(MResult,chunk) //num_threads(8)
	{
		#pragma omp for schedule(dynamic,chunk)
			for(i=startRow; i<endRow; i++){//every row m1
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



void printMatrix(int row, int col, int *M){
	/* getting pos of matrix
		->> rowi
		0 | 1 | 2
		3 | 4 | 5
		columni
		rows -> 2
		columns -> 3
		rowi*columns+columni=pos
	*/
	for(int i=0; i<row; i++){
		for(int j=0; j<col; j++){
			if(j==col-1){
				printf("%i", M[i*col+j]);
			}else{
				printf("%i | ", M[i*col+j]);
			}	
		}
		printf("\n");
	}
}






void writeTime(float elapsed, int numRows){
		/*
			Wite the time on timesc++Parallel.txt file
			M -> Matrix, Mrow -> Matrix rows, Mcol -> Matrix columns
		*/
		FILE *f = fopen("timesc++Parallel.txt","a+");//write at end of file and set result, append
		fprintf(f,"%i	%.9f\n", numRows, elapsed);
		fclose(f);
}



int main(int argc, char *argv[]) {
	//if(argc != 2){
	//	printf("There should be 2 arguments!\n");
	//	exit(1);
	//}
	MPI_Init ( &argc, &argv );
	MPI_Status status;
	int p_id, p;
	MPI_Comm_rank ( MPI_COMM_WORLD, &p_id );//identifica el número de equipo que está corriendo
	MPI_Comm_size ( MPI_COMM_WORLD, &p );//identifica el total de equipos que se usarán
  	int M1row, M1col, M2row, M2col, *M1, *M2, *MResult, *Mtemp, stop;	
	constructMatrix("inputc++.txt", M1row, M1col, M2row, M2col, M1, M2);
	auto startTime=std::chrono::high_resolution_clock::now();

	//Header part
	if(p_id==0) {
		int stepPart=M1row/(p-1), sizeTemp=0;
		MResult=new int[M1row*M2col];
		
		for(int nodeWorkerId=1,startPart=0,endPart=stepPart ; nodeWorkerId<p ; nodeWorkerId++, endPart+=stepPart)
		{
			
			if(nodeWorkerId==p-1){
				endPart=M1row;
			}
			scanf("%d", &stop);
			printf("haciendo lo primero en el worker %i", nodeWorkerId);
			MPI_Send(&startPart, 1, MPI_INT, nodeWorkerId, 0, MPI_COMM_WORLD);
			MPI_Send(&endPart, 1, MPI_INT, nodeWorkerId, 0, MPI_COMM_WORLD);
			printf("lo envie a worker %i", nodeWorkerId);
			startPart=endPart;
		}
		
		for(int nodeWorkerId=1, startPart=0, endPart=stepPart; nodeWorkerId<p; nodeWorkerId++, endPart+=stepPart) 
		{
			if(nodeWorkerId==p-1){
				endPart=M1row;
			}

			Mtemp = new int[(endPart-startPart)*M2col];
			printf("recibiendo resultado de Worker %i", nodeWorkerId);
			MPI_Recv(&Mtemp, (endPart-startPart)*M2col, MPI_INT, nodeWorkerId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("recibi el resultado del worker %i", nodeWorkerId);
			//Fill part of MatrixTemp (operator rows) into MResult
			for(int numRow=startPart, numPos=0; numRow<endPart; numRow++){
				for(int numCol=0; numCol<M1col; numCol++, numPos++){
					MResult[numRow*M1row+numCol]=Mtemp[numPos];
				}
			}
			delete [] Mtemp;
			startPart=endPart;
		}

		printMatrix(M1row, M2col, MResult);
	}

	//Workers part
	else{
		int NodeHeaderId=0, startRow=0, endRow=0;
		printf("recibiendo valores del header, soy worker %i", p_id);
		MPI_Recv(&startRow, 1, MPI_INT, NodeHeaderId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&endRow, 1, MPI_INT, NodeHeaderId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("recibi los valores del header, soy worker %i", p_id);
		
		

		mulParallelRow(startRow, endRow, M2col, M1, M2, MResult);
		printf("enviando mi resultado, soy el worker %i", p_id);
		MPI_Send(&MResult, (endRow - startRow)*M2col, MPI_INT, NodeHeaderId, 0, MPI_COMM_WORLD);
		printf("envie el resultado, soy el worker %i", p_id);
	}

	auto endTime=std::chrono::high_resolution_clock::now();
	chrono::duration<float>  elapsed = endTime - startTime;
	writeTime(elapsed.count(), M2row);
	
	// opMatrix.printOperators();
	// opMatrix.printResult();
	return 0;
}
