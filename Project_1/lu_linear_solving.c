#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <float.h>
 #define NDIM       (30)        /* number of total dimensions.  Never changes */
 
 

/*************************************************************************/
/************************************************************************
   LU_substitution():
       Performs the forward (w/ the Lower) and backward (w/ the Upper)   
       substitutions using the LU-decomposed matrix A[][] of the original
       matrix A' of the linear equation:  A'.x = B.  Upon entry, A[][]   
       is the LU matrix, B[] is the source vector, and permute[] is the  
       array containing order of permutations taken to the rows of the LU
       matrix.  See LU_decompose() for further details. 		     
    								     
       Upon exit, B[] contains the solution x[], A[][] is left unchanged.
								     
************************************************************************/
void LU_substitution( double A[][NDIM], double B[], int permute[] , int n)
{
  int i, j ;
  // int n = NDIM;
  double tmpvar,tmpvar2;

  
  /* Perform the forward substitution using the LU matrix. 
   */
  for(i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the 
       B vector to match the permutation of the LU matrix. 
       Since only the rows above the currrent one matter for 
       this row, we can permute one at a time. 
    */
    tmpvar        = B[permute[i]];
    B[permute[i]] = B[    i     ];
    for( j = (i-1); j >= 0 ; j-- ) { 
      tmpvar -=  A[i][j] * B[j];
    }
    B[i] = tmpvar; 
  }
	   

  /* Perform the backward substitution using the LU matrix. 
   */
  for( i = (n-1); i >= 0; i-- ) { 
    for( j = (i+1); j < n ; j++ ) { 
      B[i] -=  A[i][j] * B[j];
    }
    B[i] /= A[i][i] ; 
  }
  

  /* End of LU_substitution() */

}





/*************************************************************************/
/*************************************************************************
   LU_decompose():
       Performs a LU decomposition of the matrix A using Crout's method      
       with partial implicit pivoting.  The exact LU decomposition of the    
       matrix can be reconstructed from the resultant row-permuted form via  
       the integer array permute[]                                            
                                                                             
       The algorithm closely follows ludcmp.c of "Numerical Recipes  
       in C" by Press et al. 1992.                                           
                                                                             
       This will be used to solve the linear system  A.x = B                 
                                                                             
       Returns (1) if a singular matrix is found,  (0) otherwise.            
*************************************************************************/
int LU_decompose( double A[][NDIM], int permute[], int n )
{

  const  double absmin = 1.e-30; /* Value used instead of 0 for singular matrices */

  static double row_norm[NDIM];
  double  absmax, maxtemp, mintemp;

  int i, j, k, max_row;
  //int n = NDIM;


  max_row = 0;

  /* Find the maximum elements per row so that we can pretend later
     we have unit-normalized each equation: */

  for( i = 0; i < n; i++ ) { 
    absmax = 0.;
    
    for( j = 0; j < n ; j++ ) { 
      
      maxtemp = fabs( A[i][j] ); 

      if( maxtemp > absmax ) { 
	absmax = maxtemp; 
      }
    }

    /* Make sure that there is at least one non-zero element in this row: */
    if( absmax == 0. ) { 
     fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return(1);
    }

    row_norm[i] = 1. / absmax ;   /* Set the row's normalization factor. */
  }


  /* The following the calculates the matrix composed of the sum 
     of the lower (L) tridagonal matrix and the upper (U) tridagonal
     matrix that, when multiplied, form the original maxtrix.  
     This is what we call the LU decomposition of the maxtrix. 
     It does this by a recursive procedure, starting from the 
     upper-left, proceding down the column, and then to the next
     column to the right.  The decomposition can be done in place 
     since element {i,j} require only those elements with {<=i,<=j} 
     which have already been computed.
     See pg. 43-46 of "Num. Rec." for a more thorough description. 
  */

  /* For each of the columns, starting from the left ... */
  for( j = 0; j < n; j++ ) {

    /* For each of the rows starting from the top.... */

    /* Calculate the Upper part of the matrix:  i < j :   */
    for( i = 0; i < j; i++ ) {
      for( k = 0; k < i; k++ ) { 
	A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    /* Calculate the Lower part of the matrix:  i <= j :   */

    for( i = j; i < n; i++ ) {

      for (k = 0; k < j; k++) { 
	A[i][j] -= A[i][k] * A[k][j];
      }

      /* Find the maximum element in the column given the implicit 
	 unit-normalization (represented by row_norm[i]) of each row: 
      */
      maxtemp = fabs(A[i][j]) * row_norm[i] ;

      if( maxtemp >= absmax ) {
	absmax = maxtemp;
	max_row = i;
      }

    }

    /* Swap the row with the largest element (of column j) with row_j.  absmax
       This is the partial pivoting procedure that ensures we don't divide
       by 0 (or a small number) when we solve the linear system.  
       Also, since the procedure starts from left-right/top-bottom, 
       the pivot values are chosen from a pool involving all the elements 
       of column_j  in rows beneath row_j.  This ensures that 
       a row  is not permuted twice, which would mess things up. 
    */
    if( max_row != j ) {

      /* Don't swap if it will send a 0 to the last diagonal position. 
	 Note that the last column cannot pivot with any other row, 
	 so this is the last chance to ensure that the last two 
	 columns have non-zero diagonal elements.
       */

      if( (j == (n-2)) && (A[j][j+1] == 0.) ) {
	max_row = j;
      }
      else { 
	for( k = 0; k < n; k++ ) { 

	  maxtemp       = A[   j   ][k] ; 
	  A[   j   ][k] = A[max_row][k] ;
	  A[max_row][k] = maxtemp; 

	}

	/* Don't forget to swap the normalization factors, too... 
	   but we don't need the jth element any longer since we 
	   only look at rows beneath j from here on out. 
	*/
	row_norm[max_row] = row_norm[j] ; 
      }
    }

    /* Set the permutation record s.t. the j^th element equals the 
       index of the row swapped with the j^th row.  Note that since 
       this is being done in successive columns, the permutation
       vector records the successive permutations and therefore
       index of permute[] also indexes the chronology of the 
       permutations.  E.g. permute[2] = {2,1} is an identity 
       permutation, which cannot happen here though. 
    */

    permute[j] = max_row;

    if( A[j][j] == 0. ) { 
      A[j][j] = absmin;
    }


  /* Normalize the columns of the Lower tridiagonal part by their respective 
     diagonal element.  This is not done in the Upper part because the 
     Lower part's diagonal elements were set to 1, which can be done w/o 
     any loss of generality.
  */
    if( j != (n-1) ) { 
      maxtemp = 1. / A[j][j]  ;
      
      for( i = (j+1) ; i < n; i++ ) {
	A[i][j] *= maxtemp;
      }
    }

  }

  return(0);

  /* End of LU_decompose() */

}


/*************************************************************************/
/*************************************************************************
   invert_matrix():
       Uses LU decomposition and back substitution to invert a matrix 
       A[][] and assigns the inverse to Ainv[][].  This routine does not 
       destroy the original matrix A[][].                                     
                                                                             
       Returns (1) if a singular matrix is found,  (0) otherwise.            
*************************************************************************/

int invert_matrix( double Am[][NDIM], double Aminv[][NDIM], int n )  
{ 
  int i,j;
  // int n = NDIM;
  int permute[NDIM]; 
  double dxm[NDIM], Amtmp[NDIM][NDIM];

  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  Amtmp[0][i] = Am[0][i]; }

  // Get the LU matrix:
  if( LU_decompose( Amtmp,  permute, n ) != 0  ) { 
    fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
    return(1);
  }

  for( i = 0; i < n; i++ ) { 
    for( j = 0 ; j < n ; j++ ) { dxm[j] = 0. ; }
    dxm[i] = 1.; 
    
    /* Solve the linear system for the i^th column of the inverse matrix: :  */
    LU_substitution( Amtmp,  dxm, permute, n );

    for( j = 0 ; j < n ; j++ ) {  Aminv[j][i] = dxm[j]; }

  }

  return(0);
}

/******************************************
*    MULTIPLICATION OF TWO SQUARE REAL    *                                     
*    MATRICES                             *
* --------------------------------------- *		                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*N                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *		                                     
* OUTPUTS:   C  MATRIX N*N PRODUCT A*B    *                                     
*                                         *
******************************************/
void MatMult(int n, double A[][NDIM],double B[][NDIM], double C[][NDIM]) {
  double SUM;
  int I,J,K;
  for (I=0; I<n; I++)                                                                  
    for (J=0; J<n; J++) {
      SUM = 0.0;                                                                
      for (K=0; K<n; K++)
       SUM += A[I][K]*B[K][J];                                               
      C[I][J]=SUM;                                                            
    }                                                                   
}

// print a square real matrix A of size n with caption s
// (n items per line).
void MatPrint(char *s, int n, double A[][NDIM]) {
    int i,j;
    printf("\n %s\n", s);
	for (i=0; i<n; i++) {
	  for (j=0; j<n; j++) 
	     printf(" %10.6f",A[i][j]);
      printf("\n");
    }
}





/******************************************
*    MULTIPLICATION OF TWO SQUARE REAL    *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*1                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *	                                     
* OUTPUTS:   C  MATRIX N*1 PRODUCT A*B    *                                     
*                                         *
******************************************/
void MatMultOneDim(int n, double A[NDIM][NDIM],double B[], double C[]) {
  double SUM=0.0;
  int I,J,K;
  for (I=0; I<n; I++)    {                                                              
    for (J=0; J<n; J++) {

       SUM += A[I][J]*B[J];                                                              
    }                
	C[I]=SUM;  
	SUM = 0.0;
  }	
}


/******************************************
*    SUBTRACTION OF TWO REAL    	  *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*1                *                                     
*            B  MATRIX N*1                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*1 = A-B    	  *                                     
*                                         *
******************************************/
void MatSubtraction(int n, double A[] ,double B[], double C[]){
  int I,J;                                                            
    for (J=0; J<n; J++) {
                                                                      
       C[J] = A[J]-B[J];  
		if(C[J]<0) C[J]=-C[J];
                                                           
    }                

}



/******************************************
*    SUBTRACTION OF TWO SQUARE REAL   	  *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*N                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*N PRODUCT A-B	  *                                     
*                                         *
******************************************/
void MatSubtraction2(int n, double A[NDIM][NDIM] , double B[NDIM][NDIM], double C[NDIM][NDIM]){
  int I,J;
  for (I=0; I<n; I++)    {                                                              
    for (J=0; J<n; J++) {

       C[I][J] = A[I][J]-B[I][J];  
                                                      
    }                

  }	
}


/******************************************
*    SUM OF ROWS OF SQUARE REAL    	  *                                     
*    MATRIX                               *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*                  N  INTEGER             *                                     
*                             		  *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*1        	  *                                     
*                                         *
******************************************/
void MatRowsSum(int n, double A[NDIM][NDIM],double C[]){
  double SUM=0.0;
  int I,J,K;
  for (I=0; I<n; I++)    {                                                              
    for (J=0; J<n; J++) {
       SUM+= A[I][J];  
                                                       
    }                
	C[I]=SUM;  
	SUM = 0.0;
  }	
}




int readmatrix(size_t rows, size_t cols,double a[NDIM][NDIM] , const char* filename)
{

    FILE *pf;
    pf = fopen (filename, "r");
    if (pf == NULL)
        return 0;

    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t j = 0; j < cols; ++j)
            fscanf(pf, "%lf", a[i] + j);
    }


    fclose (pf); 
    return 1; 
}



int areDoubleEqual(double a, double b) {
    if (fabs(a - b) <= 1000 * DBL_EPSILON * fabs(a + b)) {
        /*A and B are equal*/
        return 1;
    }
    return 0;
}



/*This will sort in descending order. If you want ascending order, interchange 1 and -1*/
int compareDoubles(const void *a, const void *b) {
    double doubleA = *(double *) a;
    double doubleB = *(double *) b;
    if(areDoubleEqual(doubleA, doubleB)) {
        /*When A and B are equal, quick sort expects 0.
        For further details check documentation*/
        return 0;
    }
    /*A is bigger*/
    if((doubleA - doubleB) >
        ((fabs(doubleA) < fabs(doubleB) ? fabs(doubleB) :
        fabs(doubleA)) * DBL_EPSILON * 1000)) {
        return -1;
    }
    /*B is bigger*/
    return 1;
}


// main program to demonstrate the uses of functions
int main() {
  double A[NDIM][NDIM], C[NDIM][NDIM], Tol=0.01, Am1[NDIM][NDIM], Am[NDIM][NDIM], Aminv[NDIM][NDIM], B[NDIM]={25.0, 75.0, 37.0, 46.0, 5.0, 93.0, -16.0, 41.0}, B1[NDIM],B2[NDIM], arr[NDIM][NDIM];
  double C1[NDIM], C2[NDIM], x[NDIM]={1.0,  2.0,  3.0,  4.0,  4.0,  3.0,  2.0,  1.0};
  int i,j, n = 8, flag = 1, P[NDIM], permute[NDIM], choice, c=0,y=1;
  char answer;
  printf(" Inversion of a square real symetric matrix by Cholesky method\n");
  printf(" (The matrix must positive def.).\n");
  printf("\nSolution of exercise 1.a.2\n");
  printf("____________________________________________\n\n\n");
  printf("                                //////////////////////////////////////////\n");
  printf("                                //      1.Read from file 1a2.txt        //\n");
  printf("                                /////////////////////////////////////////\n");   	
  printf("\n\n\nSize of the matrix is 4:   ");
  printf("\n\n\n\n2. Read drom file 1a2.txt\n");
  printf("2. Insert random data\n");
  printf("\n\n\n\n2. Read drom file 1a2.txt\n");
  readmatrix(n, n, Am1, "1a2.txt");
  printf("The MATRIX is:\n\n\n\n");
  for(i = 0; i < n; ++i)
  {
	printf("\t|");
	for(j = 0; j < n; ++j)
		printf("\t%0.0f\t ", Am1[i][j]);	
	printf("|\n");
  }
			
  for (i=0; i<n; i++) {
	for (j=0; j<n; j++) 
		A[i][j]=Am1[i][j];					  
  }

  printf("\nSolve Ax=B linear system with LU method\n");
  printf("\n\n");

 //////////////////////                 Solve Ax=B linear system with LU method                 ////////////////////////////////////////////
           
					  
  /*Give the values of b(of Ax=b)*/
  for (i=0; i<n; i++) 
	B2[i]=B[i];
		
  LU_decompose( Am1,  permute, n );
  MatPrint("Verification A ",n,Am1);
  
  LU_substitution( Am1,  B, permute, n );
  printf("\n\n");
  for (i=0; i<n; i++)
	printf("%10.6f",B[i]);
  printf("\n\n");



																												/* Ypologismos Sfalmatwn*/	
					
//i]						
  /*We calculate the subtraction δx=||x-x'|| and place the results in the C2 matrix. Then, we sort them using the qsort algorithm to find the maximum element in C2, which becomes the numerator of the fraction.*/
  MatSubtraction(n, x, B, C2);
  qsort(C2, n, sizeof(double), compareDoubles);                 


  /*We sort the elements of the solution x using the qsort algorithm, finding their absolute values, in order to locate the maximum element and use it as the denominator of the fraction*/
  for (i=0; i<n; i++) 
     x[i]=fabs(x[i]);  
    
  qsort(x, n, sizeof(double), compareDoubles);
  printf("\nThe absolute error is:\n");		
  printf("\n%10.6f / %10.6f",C2[0], x[0]);	
  printf("\n");
		

//ii]		
  /*We apply the multiplication A*x' (where x' is the current value of the solution) and place the results in the matrix C1. */
  MatMultOneDim(n,A,B,C1);					
							
  /*We perform the subtraction δx=||b-Ax'|| and place the results in the matrix C2. Then, we sort them using the qsort algorithm in order to find the maximum element in C2, which becomes the numerator of the fraction*/
  MatSubtraction(n, B2, C1, C2);
  qsort(C2, n, sizeof(double), compareDoubles);							
  printf("\nThe absolute remainder is:\n");
  printf("\n%10.6f / %10.6f\n\n",C2[0], x[0]);
  printf("\n");

}







