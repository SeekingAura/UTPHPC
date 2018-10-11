#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
using namespace std;




class matrix{
	public:
	int *M1, *M2, *MResult; /* matrices (M1,M2) */
	size_t M1row=0,M1col=0, M2row=0, M2col=0;
	//constructor
	matrix(const char *fileName1){
		FILE *f1=NULL; /* file pointers */
		//this->buildMatrix();
		//printf("constructor\n");
		f1=this->openFile(fileName1);

		this->M1=this->buildMatrix(f1, this->M1row, this->M1col);
		this->getData(f1, this->M1, this->M1row*this->M1col);

		this->M2=this->buildMatrix(f1, this->M2row, this->M2col);
		this->getData(f1, this->M2, this->M2row*this->M2col);

		this->MResult=(int *)malloc(this->M1row*this->M2row*sizeof(int));
		
		fclose(f1);
		//f2=this->openFile(fileName2);


		// fscanf(f,"%i",this->M2row); /* %zu zx is size_t */
		// fscanf(f,"%i",this->M2col);	

	}
	//destructor
	~matrix(){
		free(this->M1);
		free(this->M2);
		free(this->MResult);
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
		int *M = (int *)malloc(rows*columns*sizeof(int)); /* reserved memory */
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
	void mulSecuential(){
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
				this->MResult[i*this->M1col+j] = data;
			}
		}
	}

	void mulParallel(size_t pos){
		size_t result=0, filaPos=pos/this->M1row, columPos=pos%this->M1row;

		for(size_t i=0;i<this->M1row; i++){
			result=result+this->M1[filaPos*this->M1row+i]*this->M2[i*this->M2col+columPos];
		}
		this->MResult[this->M1row*filaPos+columPos]=result;
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



};

void writeTime(double elapsed, size_t len){
		/*
			Wite the result on output.txt file
			M -> Matrix, Mrow -> Matrix rows, Mcol -> Matrix columns
		*/
		FILE *f = fopen("timesc++Secuential.txt","a+");//write at end of file and set result, append
		//float value=;
		fprintf(f,"%i	%.9f\n", len, elapsed);
		fclose(f);
}



int main(int argc, char const *argv[]) {
	if(argc != 2){
		printf("There should be 2 arguments!\n");
		exit(1);
	}
	matrix opMatrix(argv[1]);



	double startTime = omp_get_wtime();
	opMatrix.mulSecuential();
	double endTime = omp_get_wtime();
	double elapsed = endTime - startTime;
	writeTime(elapsed, opMatrix.M1row);
	
	// opMatrix.printOperators();
	// opMatrix.printResult();
	return 0;
}