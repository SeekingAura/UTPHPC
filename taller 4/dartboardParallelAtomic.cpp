#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <omp.h>
#include <chrono>


using namespace std;


double calcPiDartboardSecuential(long n){
	long hits=0;
	const double factor = 1.0 / RAND_MAX; //granted max value of rands 1.0
	
	double x, y;
	for(long k=0; k<n; ++k){
		/* find two numbers within 0..1*/
		x=rand()*factor;

		y=rand()*factor;        
		if(x*x + y*y < 1.0){
			++hits;
		}
	}
	return 4.0*hits/n;
}

double calcPiDartboardParallel(long n){
	long hits=0, k, chunk=n/8;
	const double factor = 1.0 / RAND_MAX;
	
	double x, y;
	int threadmasterId = 0, numThreads = 0;
	#pragma omp parallel private(x,y,k,threadmasterId,numThreads) shared(hits,chunk) //num_threads(8)
	{
		threadmasterId = omp_get_thread_num();
		if(threadmasterId == 0){
			numThreads = omp_get_num_threads();
			// printf("Thread master! with id %d\n",threadmasterId);
			// printf("Invoked threads number %d",numThreads);
		}
        // #pragma omp for nowait
        for(k=0; k<n; ++k){
            /* find two numbers within 0..1*/
            x=rand()*factor;
            y=rand()*factor;        
            if(x*x + y*y < 1.0){
                #pragma omp atomic
                    ++hits;
            }
        }
		
	}
	return 4.0*hits/n;
}



void writeTime(float elapsed, long len){
	/*
		Wite the result on output.txt file
		M -> Matrix, Mrow -> Matrix rows, Mcol -> Matrix columns
	*/
	FILE *f = fopen("timesc++ParallelAtomic.txt","a+");//write at end of file and set result, append
	//float value=;
	fprintf(f,"%ld	%.9f\n", len, elapsed);
	fclose(f);
}

int main (int argc, char const *argv[]){
	//printf("Enter the no of tosses: ");
	//cin>>n;
	if(argc != 2){
		printf("There should be 2 arguments!\n");
		exit(1);
	}
	srand((int)clock());
	long n=stol(argv[1], nullptr);
	auto startTime=chrono::high_resolution_clock::now();
	double pi_approx = calcPiDartboardParallel(n);
	auto endTime=chrono::high_resolution_clock::now();
	chrono::duration<float>  elapsed = endTime - startTime;
	writeTime(elapsed.count(), n);
	// printf("Approximation of pi after %ld tosses %.20f %.9f %% \n", n, pi_approx, fabs(M_PI-pi_approx)*100/M_PI);
	return 0;
}