#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <pthread.h>
/* --
 * Do nsweeps sweeps of Jacobi iteration on a 1D Poisson problem
 * 
 *    -u'' = f
 *
 * discretized by n+1 equally spaced mesh points on [0,1].
 * u is subject to Dirichlet boundary conditions specified in
 * the u[0] and u[n] entries of the initial vector.
 */
pthread_mutex_t mutex;
struct threadData{
	double *uWrite, *uRead, *f, h2;
    int start, end;
};

void *calcMatrixU(void *arg){
    //utmp[i] = (u[i-1] + u[i+1] + h2*f[i])/2;
	struct threadData *dataIn;
	dataIn = (struct threadData *) arg;
    printf("calculando MAtrixU\n");
    for(int i=dataIn->start; dataIn->start<dataIn->end;i++){
        pthread_mutex_lock(&mutex);
        dataIn->uWrite[i] = (dataIn->uRead[i-1] + dataIn->uRead[i+1] + dataIn->h2*dataIn->f[i])/2;
        pthread_mutex_unlock(&mutex);
    }
	//pthread_exit(NULL);

}

void jacobi(int nsweeps, int n, double* u, double* f){
    double h  = 1.0 / n;
    double h2 = h*h;
    double* utmp;
    utmp = new double[n+1];

    /* Fill boundary conditions into utmp */
    utmp[0] = u[0];
    utmp[n] = u[n];
    int numThreads=8;
    int chunk=n/numThreads;
    printf("formando hilos\n");
    pthread_t threads[numThreads];
    printf("formando data de hilos\n");
    struct threadData data[numThreads];
    pthread_mutex_init(&mutex, NULL);
    for (int sweep = 0; sweep < nsweeps; sweep += 2) {  
        /* Old data in u; new data in utmp */
        for (int i = 0, workStart=1, workEnd=chunk+1; i < numThreads; i++, workStart+=chunk, workEnd+=chunk){
            printf("creando el hilo num %i\n", i);
            data[i].uWrite=utmp;
            data[i].uRead=u;
            data[i].h2=h2;
            data[i].f=f;
            data[i].start=workStart;
            data[i].end=workEnd;
            if(i==numThreads-1 &&  workEnd!=n){
                data[i].end=n;
            }
            // utmp[i] = (u[i-1] + u[i+1] + h2*f[i])/2;
            pthread_create(&(threads[i]), NULL, calcMatrixU, (void *)&data[i]);
        }

        for(int i=0; i<numThreads; i++){
            printf("join del hilo num %i\n", i);
		    pthread_join(threads[i], NULL);
	    }
        
        
            /* Old data in utmp; new data in u */
        for (int i = 0, workStart=1, workEnd=chunk+1; i < numThreads; i++, workStart+=chunk, workEnd+=chunk){
            data[i].uWrite=u;
            data[i].uRead=utmp;
            data[i].h2=h2;
            data[i].f=f;
            data[i].start=workStart;
            data[i].end=workEnd;
            if(i==numThreads-1 &&  workEnd!=n){
                data[i].end=n;
            }
            // u[i] = (utmp[i-1] + utmp[i+1] + h2*f[i])/2;
            pthread_create(&(threads[i]), NULL, calcMatrixU, (void *)&data[i]);
        }

        for(size_t i=0; i<numThreads; i++){
		    pthread_join(threads[i], NULL);
	    }
        
        
    }
    pthread_mutex_destroy(&mutex);

    delete [] utmp;
}


void write_solution(int n, double* u, const char* fname){
    double h = 1.0 / n;
    FILE* fp = fopen(fname, "w+");
    for (int i = 0; i <= n; ++i){
        fprintf(fp, "%g %g\n", i*h, u[i]);
    }
    fclose(fp);
}

void writeTime(float elapsed, int valueK, int len){
	FILE *f = fopen("timesc++Thread.txt","a+");//write at end of file and set result, append
	fprintf(f,"%ld	%ld	%.9f\n", valueK, len, elapsed);
	fclose(f);
}

int main(int argc, char** argv){
    if(argc < 2){
		printf("There should be 3 or 2 arguments! (outputfile is optional)-> nameexecute.exe $n $steps $outputfile\n");
		exit(1);
	}
    int n, nsteps;
    double* u;
    double* f;
    double h;
    char* fname;

    /* Process arguments */
    n      = (argc > 1) ? atoi(argv[1]) : 100;
    nsteps = (argc > 2) ? atoi(argv[2]) : 100;
    fname  = (argc > 3) ? argv[3] : NULL;
    h      = 1.0/n;

    /* Allocate and initialize arrays */
    u = new double[n+1];
    f = new double[n+1];
    memset(u, 0, (n+1) * sizeof(double));
    for (int i = 0; i <= n; ++i){
        f[i] = i * h;
    }

    /* Run the solver */
    printf("empezó\n");
    auto startTime=std::chrono::high_resolution_clock::now();
    jacobi(nsteps, n, u, f);
    auto endTime=std::chrono::high_resolution_clock::now();
    std::chrono::duration<float>  elapsed = endTime - startTime;
    printf("termino\n");
	writeTime(elapsed.count(), n, nsteps);
    //get_time(&tend);

    /* Run the solver */    
    /*printf("n: %d\n"
           "nsteps: %d\n"
           "Elapsed time: %g s\n", 
           n, nsteps, timespec_diff(tstart, tend));
*/
    /* Write the results */
    if (fname){
        write_solution(n, u, fname);
    }

    delete [] f;
    delete [] u;
    return 0;
}
