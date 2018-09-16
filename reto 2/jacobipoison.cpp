#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <chrono>
#include <omp.h>
using namespace std;



void r8vec_print ( long n, double *a){

    /*
        ***************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT prints an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, double A[N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
    */

    printf("\n");
    printf("Solución");
    printf("\n");
    printf("\n");
    for (long i = 0; i < n; i++ ){
        printf("  %ld", i);
        printf(": %.9f", a[i]);  
        printf("\n");
    }

    return;
}

double* linspace(double startValue, double endValue, long n){
    double factor=endValue-startValue, *result;
    result= new double[n];
    if(factor<0){
        printf("ERROR factor negative");
    }
    factor=factor/(n-1);
    //long chunk=n/8;
    //#pragma omp parallel private(i) shared(chunk, result, factor, startValue) //num_threads(8)
	//{
        //#pragma omp for nowait
        for(long i=0;i<n;i++){
            result[i]=startValue+factor*i;
        }
    //}
    return result;

}


double* force(long sizeVector){
    double *vectorResult, factor;
    vectorResult=new double[sizeVector];
    // long chunk=sizeVector/8;
    //#pragma omp parallel private(i, factor) shared(chunk, vectorResult) //num_threads(8)
    //{
        //#pragma omp for schedule(dynamic, chunk)
            for(long i=0; i<sizeVector; i++){
                factor=((double) i)/((double)(sizeVector-1));//equivalent to linespace 0 to 1 in n
                vectorResult[i]=(-factor)*(factor+3)*(exp(factor));
            }
    //}
    return vectorResult;
    
}

double* exact(long sizeVector){
    double *vectorResult, factor;
    vectorResult=new double[sizeVector];
    //long chunk=sizeVector/8;
    #pragma omp parallel private(i, factor) shared(chunk, vectorResult, sizeVector) //num_threads(8)
    {
        #pragma omp for schedule(dynamic, chunk)
        for(long i=0; i<sizeVector; i++){

            factor=((double) i)/((double)(sizeVector-1));
            vectorResult[i]=(factor)*(factor-1)*(exp(factor));
        }
    }
    return vectorResult;
    
}

long potenciaLong(long base, long exponente){
    long result=base;
    if(exponente==0){
        return 1;
    }
    //#pragma omp parallel for reduction(*:result)
    for(long i=1; i<exponente; i++){
        result*=base;
    }
    return result;
}

double potenciaDouble(double base, long exponente){
    double result=base;
    if(exponente==0){
        return 1;
    }
    //#pragma omp parallel for reduction(*:result)
    for(long i=1; i<exponente; i++){
        result*=base;
    }
    return result;
}

double *matrixMakeA ( long m, long n, double hk){

    /*
        ***************************************************************************80
        //
        //  Purpose:
        //
        //    DIF2 returns the DIF2 matrix.
        //
        //  Example:
        //
        //    N = 5
        //
        //    2 -1  .  .  .
        //   -1  2 -1  .  .
        //    . -1  2 -1  .
        //    .  . -1  2 -1
        //    .  .  . -1  2
        //
        //  Properties:
        //
        //    A is banded, with bandwidth 3.
        //
        //    A is tridiagonal.
        //
        //    Because A is tridiagonal, it has property A (bipartite).
        //
        //    A is a special case of the TRIS or tridiagonal scalar matrix.
        //
        //    A is integral, therefore det ( A ) is integral, and 
        //    det ( A ) * inverse ( A ) is integral.
        //
        //    A is Toeplitz: constant along diagonals.
        //
        //    A is symmetric: A' = A.
        //
        //    Because A is symmetric, it is normal.
        //
        //    Because A is normal, it is diagonalizable.
        //
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //
        //    A is positive definite.
        //
        //    A is an M matrix.
        //
        //    A is weakly diagonally dominant, but not strictly diagonally dominant.
        //
        //    A has an LU factorization A = L * U, without pivoting.
        //
        //      The matrix L is lower bidiagonal with subdiagonal elements:
        //
        //        L(I+1,I) = -I/(I+1)
        //
        //      The matrix U is upper bidiagonal, with diagonal elements
        //
        //        U(I,I) = (I+1)/I
        //
        //      and superdiagonal elements which are all -1.
        //
        //    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
        //
        //      L(I,I) =    sqrt ( (I+1) / I )
        //      L(I,I-1) = -sqrt ( (I-1) / I )
        //
        //    The eigenvalues are
        //
        //      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
        //                = 4 SIN^2(I*PI/(2*N+2))
        //
        //    The corresponding eigenvector X(I) has entries
        //
        //       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
        //
        //    Simple linear systems:
        //
        //      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
        //
        //      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
        //
        //    det ( A ) = N + 1.
        //
        //    The value of the determinant can be seen by induction,
        //    and expanding the determinant across the first row:
        //
        //      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
        //                = 2 * N - (N-1)
        //                = N + 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Gregory, David Karney,
        //    A Collection of Matrices for Testing Computational Algorithms,
        //    Wiley, 1969,
        //    ISBN: 0882756494,
        //    LC: QA263.68
        //
        //    Morris Newman, John Todd,
        //    Example A8,
        //    The evaluation of matrix inversion programs,
        //    Journal of the Society for Industrial and Applied Mathematics,
        //    Volume 6, Number 4, pages 466-476, 1958.
        //
        //    John Todd,
        //    Basic Numerical Mathematics,
        //    Volume 2: Numerical Algebra,
        //    Birkhauser, 1980,
        //    ISBN: 0817608117,
        //    LC: QA297.T58.
        //
        //    Joan Westlake,
        //    A Handbook of Numerical Matrix Inversion and Solution of 
        //    Linear Equations,
        //    John Wiley, 1968,
        //    ISBN13: 978-0471936756,
        //    LC: QA263.W47.
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double DIF2[M*N], the matrix.
        //
    */
    double *a;
    a=new double[m*n];
    double divValue=potenciaDouble(hk, 2);
    //long chunk=n/8;


    //#pragma omp parallel private(i, j) shared(chunk, a, divValue, m, n) //num_threads(8)
    //{
        //#pragma omp for schedule(dynamic, chunk)
        for (long j = 0; j < n; j++ ){
            for (long i = 0; i < m; i++ ){
                if ( j == i - 1 ){
                    a[i+j*m] = ((double)-1.0)/divValue;
                }
                else if ( j == i ){
                    a[i+j*m] = ((double)2.0)/divValue;
                }
                else if ( j == i + 1 ){
                    a[i+j*m] = ((double)-1.0)/divValue;
                }
                else{
                    a[i+j*m] = (double)0.0;
                }
            }
        }
    //}
    a[0]=1.0;
    a[m*n-1]=1.0;
    
    return a;
}

double *matrixMakeA ( long m, long n){

    /*
        ***************************************************************************80
        //
        //  Purpose:
        //
        //    DIF2 returns the DIF2 matrix.
        //
        //  Example:
        //
        //    N = 5
        //
        //    2 -1  .  .  .
        //   -1  2 -1  .  .
        //    . -1  2 -1  .
        //    .  . -1  2 -1
        //    .  .  . -1  2
        //
        //  Properties:
        //
        //    A is banded, with bandwidth 3.
        //
        //    A is tridiagonal.
        //
        //    Because A is tridiagonal, it has property A (bipartite).
        //
        //    A is a special case of the TRIS or tridiagonal scalar matrix.
        //
        //    A is integral, therefore det ( A ) is integral, and 
        //    det ( A ) * inverse ( A ) is integral.
        //
        //    A is Toeplitz: constant along diagonals.
        //
        //    A is symmetric: A' = A.
        //
        //    Because A is symmetric, it is normal.
        //
        //    Because A is normal, it is diagonalizable.
        //
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //
        //    A is positive definite.
        //
        //    A is an M matrix.
        //
        //    A is weakly diagonally dominant, but not strictly diagonally dominant.
        //
        //    A has an LU factorization A = L * U, without pivoting.
        //
        //      The matrix L is lower bidiagonal with subdiagonal elements:
        //
        //        L(I+1,I) = -I/(I+1)
        //
        //      The matrix U is upper bidiagonal, with diagonal elements
        //
        //        U(I,I) = (I+1)/I
        //
        //      and superdiagonal elements which are all -1.
        //
        //    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
        //
        //      L(I,I) =    sqrt ( (I+1) / I )
        //      L(I,I-1) = -sqrt ( (I-1) / I )
        //
        //    The eigenvalues are
        //
        //      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
        //                = 4 SIN^2(I*PI/(2*N+2))
        //
        //    The corresponding eigenvector X(I) has entries
        //
        //       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
        //
        //    Simple linear systems:
        //
        //      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
        //
        //      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
        //
        //    det ( A ) = N + 1.
        //
        //    The value of the determinant can be seen by induction,
        //    and expanding the determinant across the first row:
        //
        //      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
        //                = 2 * N - (N-1)
        //                = N + 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Gregory, David Karney,
        //    A Collection of Matrices for Testing Computational Algorithms,
        //    Wiley, 1969,
        //    ISBN: 0882756494,
        //    LC: QA263.68
        //
        //    Morris Newman, John Todd,
        //    Example A8,
        //    The evaluation of matrix inversion programs,
        //    Journal of the Society for Industrial and Applied Mathematics,
        //    Volume 6, Number 4, pages 466-476, 1958.
        //
        //    John Todd,
        //    Basic Numerical Mathematics,
        //    Volume 2: Numerical Algebra,
        //    Birkhauser, 1980,
        //    ISBN: 0817608117,
        //    LC: QA297.T58.
        //
        //    Joan Westlake,
        //    A Handbook of Numerical Matrix Inversion and Solution of 
        //    Linear Equations,
        //    John Wiley, 1968,
        //    ISBN13: 978-0471936756,
        //    LC: QA263.W47.
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double DIF2[M*N], the matrix.
        //
    */

    double *a;
    a=new double[m*n];
    
    //long chunk=n/8; 
    //#pragma omp parallel private(i, j) shared(chunk, a, n, m) //num_threads(8)
    //{
    
        //#pragma omp for schedule(dynamic, chunk)
        for (long j = 0; j < n; j++ ){
            
            for (long i = 0; i < m; i++ ){
                if ( j == i - 1 ){
                    a[i+j*m] = (double)-1.0;
                    //printf("value a[i+j*m]%.9f", a[i+j*m]);
                }
                else if ( j == i ){
                    a[i+j*m] = (double)2.0;
                }
                else if ( j == i + 1 ){
                    a[i+j*m] = (double)-1.0;
                }
                else{
                    a[i+j*m] = 0.0;
                }
            }
        }
    //}

    return a;
}

double *r8mat_mv_new ( long m, long n, double *a, double *x ){

    /*
        ***************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MV_NEW multiplies a matrix times a vector.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8MAT_MV_NEW[M], the product A*X.
        //
    */

    double *y;
    y=new double[m];
    //long chunk=m/8;
    //#pragma omp parallel private(i, j) shared(chunk, a, x, y, m, n) //num_threads(8)
    //{
    
        //#pragma omp for schedule(dynamic, chunk)
        
        for (long i = 0; i < m; i++ ){
            y[i] = 0.0;
            for (long j = 0; j < n; j++ ){
                y[i] = y[i] + a[i+j*m] * x[j];
            }
        }
    //}
    return y;
}

double * zeros(long n){
    double *x;
    x=new double[n];
    //long chunk=n/8;
    //#pragma omp parallel private(i) shared(chunk, x, n) //num_threads(8)
    //{
    
        //#pragma omp for nowait
        for (long i = 0; i < n; i++ ){
            x[i] = 0.0;
        }
    //}
    return x;
}

double r8mat_residual_norm ( long m, long n, double *a, double *x, double *b){

    /*
        ***************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_RESIDUAL_NORM returns the norm of A*x-b.
        //
        //  Discussion:
        //
        //    A is an MxN R8MAT, a matrix of R8's.
        //
        //    X is an N R8VEC, and B is an M R8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the right hand side vector.
        //
        //    Output, double R8MAT_RESIDUAL_NORM, the norm of A*x-b.
        //
    */

    double *r;
    r=new double[m];
    double r_norm;

    //long chunk=m/8;
    //#pragma omp parallel private(i, j) shared(chunk, r, a, b, x, m, n) //num_threads(8)
    //{
    
        //#pragma omp for nowait
        for (long i = 0; i < m; i++ ){
            r[i] = - b[i];
            //#pragma omp parallel for reduction(+:r[i])
            for (long j = 0; j < n; j++ ){
                r[i] +=a[i+j*m] * x[j];
            }
        }
    //}

    r_norm = 0.0;
    //#pragma omp parallel for reduction(+:r_norm)
    for (long i = 0; i < m; i++ ){
        r_norm+=r[i] * r[i];
    }
    r_norm = sqrt ( r_norm );

    delete [] r;

    return r_norm;
}

double r8vec_diff_norm_squared ( long n, double *a, double *b ){

    /*
        ***************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The square of the L2 norm of the difference of A and B is:
        //
        //      R8VEC_DIFF_NORM_SQUARED = sum ( 1 <= I <= N ) ( A[I] - B[I] )^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], B[N], the vectors.
        //
        //    Output, double R8VEC_DIFF_NORM_SQUARED, the square of the L2 norm of A - B.
        //
    */

    double value= 0.0;
    //No paralelizable
    //#pragma omp for nowait
    for (long i = 0; i < n; i++ ){
        value += ( a[i] - b[i] ) * ( a[i] - b[i] );
    }

    return value;
}

double *jacobi1 ( long n, double *a, double *b, double *x){

    /*
        ***************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI1 carries out one step of the Jacobi iteration.
        //
        //  Discussion:
        //
        //    The linear system A*x=b is to be solved.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N,N], the matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Input, double X[N], the current solution estimate.
        //
        //    Output, double JACOBI1[N], the solution estimate updated by
        //    one step of the Jacobi iteration.
        //
    */
    double *x_new;
    x_new=new double[n];

    // long chunk=n/8;
    //#pragma omp parallel private(i, j) shared(chunk, n, x_new, b, a, x) //num_threads(8)
    //{
    
     //   #pragma omp for schedule(dynamic, chunk)
        for (long i = 0; i < n; i++ ){
            x_new[i] = b[i];
            //#pragma omp for nowait
            for (long j = 0; j < n; j++ ){
                if ( j != i ){
                    x_new[i] = x_new[i] - a[i+j*n] * x[j];
                }
            }
            x_new[i] = x_new[i] / a[i+i*n];
        }
    //}
  return x_new;
}



    


void jacobi_poisson_1d (long k, long iterations){
    
    /*
    %
    % Set boundaries.
    %
    */
    double a = 0.0;
    double b = 1.0;
    /*
    %
    % Set boundary conditions.
    %
    */
    double ua = 0.0;
    double ub = 0.0;
    /*
    %
    % Get NK.
    %(long)val 2^k+1
    */
    printf("haciendo potenciaLong\n");
    long nk = potenciaLong(2, k) + 1;
    
    /*%
    % Set XK.
    %*/
    // double *xk = linspace ( a, b, nk );
    /*%//xk is now inside on force and exact value
    % Get HK HK is w.
    %*/
    double w = ( b - a ) / ( nk - 1 );
    /*%
    % Set FK.
    %*/
    //fk=force(nk);
    //fk[0]=ua;
    //fk[nk-1]=ub;
    /*%
    % Set the -1/2/-1 entries of A.
    %
    % In order that the operator A approximation the Poisson operator,
    % and in order that we can compare linear systems for successive grids,
    % we should NOT multiply through by hk^2.
    %
    % Though it is tempting to try to "normalize" the matrix A, the
    % unlongended result is to scale our right hand side a multiplicative
    % factor of hk^2, which means that we make it easier and easier to
    % satisfy the RMS residual tolerance, as NK increases, with solution
    % vectors that are actually worse and worse.
    %*/
    printf("Obteniendo matrix A\n");
    double *A = matrixMakeA(nk, nk);
    
    // double *A = matrixMakeA(nk, nk, w);
    /*
    %
    % Just because we can, ask MATLAB to get the exact solution of the linear system
    % directly.
    %
    */
    //udk = A \ fk; JODIDO
    /*
    %
    % Sample the solution to the continuous problem.
    %
    */
    printf("Obteniendo exact\n");
    double *uek = exact(nk);
    
    /*
    %
    % Use Jacobi iteration to solve the linear system to the given tolerance.
    %
    */
    printf("Obteniendo B\n");
    double *B = r8mat_mv_new ( nk, nk, A, uek );
    
    //
    //  Set the initial estimate for the solution.
    //
    printf("llenando de ceros\n");
    double *ujk = zeros (nk);
    // double tol = 0.000001;
    printf("value iterations %ld ", iterations);
    double *m_plot;
    m_plot=new double[iterations+1];
    double *r_plot;
    r_plot=new double[iterations+1];
    double *s_plot;
    s_plot=new double[iterations+1];
    double *x_plot;
    x_plot=new double[nk];//nk*iterations+1
    if(x_plot==NULL){
        printf("OJOO\n");
    }
    printf("la primera normal\n");
    r_plot[0] = r8mat_residual_norm ( nk, nk, A, ujk, B );
    m_plot[0] = 1.0;
    //
    //  Initialize plot arrays.
    //
    printf("estableciendo x_plot\n");
    // long chunk=nk/8;
    //#pragma omp parallel private(i) shared(chunk, nk, x_plot, ujk) //num_threads(8)
    //{
        //#pragma omp for nowait
        for (long i = 0; i < nk; i++ ){
            x_plot[i] = ujk[i];
        }
    //}
    
    printf("estableciendo s_plot\n");
    // chunk=iterations/8;
    
    //#pragma omp parallel private(i) shared(chunk, iterations, s_plot) //num_threads(8)
    //{
    
        //#pragma omp for nowait
        for (long i = 0; i <= iterations; i++ ){
            s_plot[i] = ( double ) i;
        }
    //}
    //
    //  Carry out the iteration.
    //
    double *ujk_new;
    //double w = 0.5;
    printf("haciendo el jacobi\n");
    for (long it = 1; it <= iterations; it++ ){
        //printf("hago jacobi1\n");
        ujk_new = jacobi1 ( nk, A, B, ujk );
        //printf("encuentro residual norm\n");


        //#pragma omp sections
        //{
           // {
                r_plot[it] = r8mat_residual_norm ( nk, nk, A, ujk_new, B );
            //}
            //
            //  Compute the average point motion.
            //
            //printf("encuentro diff norm squared\n");
            //#pragma omp section
            //{
                m_plot[it] = r8vec_diff_norm_squared ( nk, ujk, ujk_new ) / ( double ) nk;
            //}
            //
            //  Update the solution
            //
            //printf("actualizando ujk norm\n");
            
            //#pragma omp for nowait
        //}
        for (long i = 0; i < nk; i++ ){
            ujk[i] = ( 1.0 - w ) * ujk[i] + w * ujk_new[i];
        }
                
            


        //  r8vec_copy ( n, x_new, x );
        //printf("actualizando x_plot norm\n");
        //chunk=nk/8;
        //#pragma omp parallel private(i) shared(chunk, nk, x_plot, ujk) //num_threads(8)
        //{
        
            //#pragma omp for nowait
            for (long i = 0; i < nk; i++ ){
                    x_plot[i] = ujk[i];
            }
        //}

        delete [] ujk_new;
    }

    r8vec_print ( nk, ujk);

   
}


void writeTime(float elapsed, long valueK, long len){
	/*
		Wite the result on output.txt file
		M -> Matrix, Mrow -> Matrix rows, Mcol -> Matrix columns
	*/
	FILE *f = fopen("timesc++Parallel.txt","a+");//write at end of file and set result, append
	//float value=;
	fprintf(f,"%ld	%ld	%.9f\n", valueK, len, elapsed);
	fclose(f);
}

int main(int argc, char const *argv[]) {
    if(argc != 3){
		printf("There should be 2 arguments!\n");
		exit(1);
	}
    long k=stol(argv[1], nullptr);
    long iterations=stol(argv[2], nullptr);
    auto startTime=chrono::high_resolution_clock::now();
    jacobi_poisson_1d(k, iterations);
    auto endTime=chrono::high_resolution_clock::now();
    chrono::duration<float>  elapsed = endTime - startTime;
	writeTime(elapsed.count(), k, iterations);

}