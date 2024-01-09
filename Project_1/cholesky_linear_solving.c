// /* A.M.:1115201800040 */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define  SIZE 25
#define NDIM       (10)        /* number of total dimensions.  Never changes */



/**************************************************************************************
*  Solve Ax=B linear system with cholesky method					     			                  *
*  The matrix must be positive definite.                         			    			      *	 
* --------------------------------------------------------------                      *
* 					                                                                        	*
*             																		           				                  *
*             													                                        	    *
* ------------------------------------------------------------------------------------* 
* SAMPLE RUN:                                                                       	*
*Size of the matrix is 4:   



2. Read drom file 2a1.txt
The MATRIX is:



	|	5	 	7	 	6	 	5	 |
	|	7	 	10	 	8	 	7	 |
	|	6	 	8	 	10	 	9	 |
	|	5	 	7	 	9	 	10	 |

Solve Ax=B linear system with cholesky method





These are the values of b(of Ax=b)
B[0]=23.000000
B[1]=32.000000
B[2]=33.000000
B[3]=31.000000
This is L matrix:
	     2.236	     0.000	     0.000	     0.000
	     3.130	     0.447	     0.000	     0.000
	     2.683	    -0.894	     1.414	     0.000
	     2.236	     0.000	     2.121	     0.707


This is U matrix:
	     2.236	     3.130	     2.683	     2.236
	     0.000	     0.447	    -0.894	     0.000
	     0.000	     0.000	     1.414	     2.121
	     0.000	     0.000	     0.000	     0.707


lyseis

1.000000
1.000000
1.000000
1.000000

The absolute error is:

  0.000000 /   1.000000

The absolute remainder is:

  0.000000 /   1.000000
  ****************************************************************************************/



typedef double MAT[SIZE][SIZE], VEC[SIZE];
double y2[SIZE];
void choldc1(int,MAT,VEC); 

/* -----------------------------------------------
        Cholesky decomposition.

        input    n  size of matrix
        input    A  Symmetric positive def. matrix
        output   a  lower deomposed matrix
        uses        choldc1(int,MAT,VEC)
   ----------------------------------------------- */
void choldc(int n,MAT A, MAT a) {
   int i,j;
   VEC p;
   for (i = 0; i < n; i++) 
	 for (j = 0; j < n; j++) 
	   a[i][j] = A[i][j];
   choldc1(n, a, p);
   for (i = 0; i < n; i++) {
      a[i][i] = p[i];
      for (j = i + 1; j < n; j++)
         a[i][j] = 0;
   }
}
	

 
/* -----------------------------------------------------
         Inverse of Cholesky decomposition.

         input    n  size of matrix
         input    A  Symmetric positive def. matrix
         output   a  inverse of lower deomposed matrix
         uses        choldc1(int,MAT,VEC)         
   ----------------------------------------------------- */
void choldcsl(int n, MAT A, MAT a) {
  int i,j,k; double sum;
  VEC p;
      for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++) 
      a[i][j] = A[i][j];
      choldc1(n, a, p);
      for (i = 0; i < n; i++) {
        a[i][i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
          sum = 0;
          for (k = i; k < j; k++) {
            sum -= a[j][k] * a[k][i];
      }
          a[j][i] = sum / p[j];
    }
  }
}
 
/* -----------------------------------------------------------------------------
        Computation of Determinant of the matrix using Cholesky decomposition

        input    n  size of matrix
        input    a  Symmetric positive def. matrix
        return      det(a)
        uses        choldc(int,MAT,MAT)
   ------------------------------------------------------------------------------ */
double choldet(int n, MAT a) {
   MAT c; double d=1; int i;
       choldc(n,a,c);
       for (i = 0; i < n; i++)  d *= c[i][i];
       return d * d;
}
 
/* ---------------------------------------------------
        Matrix inverse using Cholesky decomposition

        input    n  size of matrix
        input	  A  Symmetric positive def. matrix
        output   a  inverse of A
        uses        choldc1(MAT, VEC)
   --------------------------------------------------- */
void cholsl(int n, MAT A, MAT a) {
  int i,j,k;
      choldcsl(n,A,a);
      for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
          a[i][j] = 0.0;
    }
  }
      for (i = 0; i < n; i++) {
        a[i][i] *= a[i][i];
        for (k = i + 1; k < n; k++) {
          a[i][i] += a[k][i] * a[k][i];
    }
        for (j = i + 1; j < n; j++) {
          for (k = j; k < n; k++) {
            a[i][j] += a[k][i] * a[k][j];
      }
    }
  }
      for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
          a[i][j] = a[j][i];
    }
  }
}

/* ----------------------------------------------------
        main method for Cholesky decomposition.

        input         n  size of matrix
        input/output  a  Symmetric positive def. matrix
        output        p  vector of resulting diag of a
   ----------------------------------------------------- */
void choldc1(int n, MAT a, VEC p) {
  int i,j,k;
  double sum;

  for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
          sum = a[i][j];
          for (k = i - 1; k >= 0; k--) {
            sum -= a[i][k] * a[j][k];
      }
          if (i == j) {
            if (sum <= 0) {
              printf(" a is not positive definite!\n");
			}	
			// if(sum<0)
                // p[i] = sqrt(-sum);
			else 
				p[i] = sqrt(sum);
      }
          else {
            a[j][i] = sum / p[i];
      }
    }
  }
}








int func(double A[NDIM][NDIM],int n,double B[SIZE], double x[SIZE] )
{
    int i,j,k;
    double pivot = 0.0;
    double factor = 0.0;
    double sum = 0.0;
   // n=4;

    double a[SIZE][SIZE];

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            a[i+1][j+1]=A[i][j];
        }    
    }

    double b[SIZE];
    
    for(i=0;i<n;i++)
    {
            b[i+1]=B[i];
    
    }


    for(i=1;i<=n;i++)
    {
        a[i][n+1]=b[i];
    }


for(k=1;k<=n-1;k++)
{

	    if(a[k][k]==0.0)
	    {

		    printf("error");

	    }
	    else
	    {
		    pivot = a[k][k];
		    for(j=k;j<=n+1;j++)
			    a[k][j]= a[k][j]/pivot; 

		    for(i=k+1;i<=n;i++)
		    {
			    factor = a[i][k];

			    for(j = k;j<=n+1;j++)
			    {
				    a[i][j] = a[i][j] - factor * a[k][j];
			    }
		    } 
	    } 

if(a[n][n]==0)

printf("error");

else
{

    x[n] = a[n][n+1]/a[n][n];

    for(i=n-1;i>=1;i--)
    {

        sum = 0.0;

        for(j=i+1;j<=n;j++)

        sum = sum + a[i][j] * x[j];

        x[i]= ( a[i][n+1]- sum )/a[i][i];
	    
    }

}

} 

    for(i=1;i<=n;i++)
    {
	    y2[i]=x[i];
    }

}


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
void LU_substitution( double A[][NDIM], double B[], int n)
{
  int i, j ;
  //int n = NDIM;
  double tmpvar,tmpvar2;
  double a[NDIM][NDIM]={0.0}, y[NDIM]={0.0},y1[NDIM]={0.0},x[SIZE]={0.0};

  /* Perform the forward substitution using the LU matrix. */
 
 
  /*Store the lower half of the A[NDIM][NDIM] matrix in a[NDIM][NDIM]*/  
    for(i = 0; i < n; i++) {
	    for( j = (i); j >= 0 ; j-- ) { 
	      a[i][j]=  A[i][j];
	    }
    }

    
  /*Solve the linear system Ly=B*/
    func(a,n,B,y);
  
  /* Perform the backward substitution using the LU matrix. */
     for(i=0;i<=n; i++)
   {
       for(j=0; j<=n;j++)
       {
	       a[i][j]=0.0;
       }
   }
  
   /*Store the upper half of the A[NDIM][NDIM] matrix in a[NDIM][NDIM] */
    for( i = (n-1); i >= 0; i-- ) { 
	    for( j = (i); j < n ; j++ ) { 
	      a[i][j]=  A[i][j];
	    } 
    }
  // printf("\n");
      for(i=0;i<n;i++)
    {
            y[i]=y2[i+1];
    }  
    
   /*Solve the linear system L^(T)x=y*/ 
     func(a,n,y,x);
     
     printf("\n\nSolution\n");
     printf("\n");
         for(i=0;i<n;i++)
    {
            B[i]=y2[i+1];
		    printf("%lf\n",B[i]);
    }

  /* End of LU_substitution() */

}




	// print a square real matrix A of size n with caption s
	// (n items per line).
    void MatPrint(char *s, int n, MAT A) {
      int i,j;
      printf("\n %s\n", s);
	  for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) 
	      printf(" %10.6f",A[i][j]);
            printf("\n");
          }
        }

	void MatPrint1(char *s, int n, double A[][NDIM]) {
          int i,j;
          printf("\n %s\n", s);
	  for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) 
	      printf(" %10.6f",A[i][j]);
            printf("\n");
          }
        }

    //Print one dimensional array
    void printArray(double arr[], int size)
    {
        int i;
        for (i = 0; i < size; i++)
            printf("%10.6lf ", arr[i]);
        printf("\n");
    }


        // check if matrix A is positive definite (return 1)
        // or not positive definite (return 0) 
        int Check_Matrix(int n, MAT A) {
          int i,j,k,result; double sum;
          result=1;
	  for (i=0; i<n; i++) {
	    for (j = i; j<n; j++) {
              sum = A[i][j];
              for (k = i - 1; k>=0; k--)
                sum -= A[i][k] * A[j][k];
              if (i == j)
                if (sum <= 0.0) result=0;
	    }
          }
	  return result;
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
void MatMult(int n, MAT A,MAT B, MAT C) {
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
void MatMultOneDim(int n, MAT A,double B[], double C[]) {
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
*    SUBTRACTION OF TWO REAL    *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*1                *                                     
*            B  MATRIX N*1                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*1 = A-B    *                                     
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
*    SUBTRACTION OF TWO SQUARE REAL    *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*N                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*N PRODUCT A-B    *                                     
*                                         *
******************************************/
void MatSubtraction2(int n, MAT A , MAT B, MAT C){
  int I,J;
  for (I=0; I<n; I++)    {                                                              
    for (J=0; J<n; J++) {

       C[I][J] = A[I][J]-B[I][J];  
                                                      
    }                

  }	
}


/******************************************
*    SUM OF ROWS OF SQUARE REAL           *                                     
*    MATRIX                               *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*                  N  INTEGER             *                                     
*                                         *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*1                *                                     
*                                         *
******************************************/
void MatRowsSum(int n, MAT A,double C[]){
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





// copy MAT A in MAT A1
void MatCopy(int n, MAT A, MAT A1) {
  int i,j;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      A1[i][j]=A[i][j];
}






int readmatrix(size_t rows, size_t cols,MAT a , const char* filename)
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






// main program to demonstrate the use of our functions
void main() {
  MAT A,Am1;    //Am1 is our matrix and A is a copy of Am1 
  MAT C, C4;
  VEC p;
  int i,j, n, r;
  double zero=0;
  double B1[SIZE]={23.0, 32.0, 33.0, 31.0},B2[SIZE], arr[NDIM][NDIM];
  double C1[SIZE], C2[SIZE], x[SIZE]={1.0, 1.0, 1.0, 1.0};
 


  int c=0;
  printf("\nSolution of exercise 2.b.1\n");
  printf("____________________________________________\n\n\n");
  printf("                                //////////////////////////////////////////\n");
  printf("                                //      1.Read from file 2a1.txt        //\n");
  printf("                                /////////////////////////////////////////\n");             	
  printf("\n\n\nSize of the matrix is 4:   ");
  printf("\n\n\n\n2. Read drom file 2a1.txt\n");
  n=4;
  readmatrix(4, 4, Am1, "2a1.txt");
  printf("The MATRIX is:\n\n\n\n");
  
  for(size_t i = 0; i < 4; ++i)
  {
     printf("\t|");
	 for(size_t j = 0; j < 4; ++j)
		printf("\t%0.0f\t ", Am1[i][j]);
	 printf("|\n");
  }
	
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) 
	   A[i][j]=Am1[i][j];					  
  }

   printf("\nSolve Ax=B linear system with cholesky method\n");
   printf("\n\n");

 ////////////////////             Solve Ax=B linear system with cholesky method                 /////////////////////////////////////////
					  
   /*Give the values of b(of Ax=b)*/
   printf("\n\n\nThese are the values of b(of Ax=b)\n");
   for (i=0; i<n; i++) {
      printf("B[%d]=%lf\n",i,B1[i]);
	}				

   for (i=0; i<n; i++) {
	  B2[i]=B1[i];
   }
				  

	/* choldc1 function finds L matrix  */
	choldc1(n, Am1, p);
				 
	/* we  store L in arr[][]*/
	for(r = 0; r < n; r++){
		for(c = 0; c < n; c++){
			if(r > c)
				arr[r][c]=Am1[r][c];     
			if(r == c)
				arr[r][c]=p[c];
			else if(r<c)
				arr[r][c]=0.0;
		}
	}
				   

	/* Then we find L^(T) by symmetry and we store it in arr[][] too*/
	for (i=0; i<n; i++)
		for (j=i+1; j<n; j++)
			arr[i][j]=arr[j][i];
					   
	/* Let's print L and L^T matrixes respectively*/
	printf("This is L matrix:\n");
	for(i = 0; i < n; i++) {
		for( j = 0; j < n ; j++ ) { 
			if(j>i)
				printf("\t%10.3f",zero);
			else
				printf("\t%10.3f",arr[i][j]);
		}
			printf("\n");
	}
					
	printf("\n\nThis is U matrix:\n");
	for (i=0; i<n; i++){
	   for (j=0; j<n; j++){
		   if(j<i)
			   printf("\t%10.3f",zero);
			else
			   printf("\t%10.3f",arr[i][j]);
		}
						  printf("\n");
	}
					
	LU_substitution(arr, B1, n);
									
					
															/*Calculation of Errors.*/	
					
//i]						
    /*We calculate the subtraction δx=||x-x'|| and place the results in the C2 matrix. Then, we sort them using the qsort algorithm to find the maximum element in C2, which becomes the numerator of the fraction.*/
	MatSubtraction(n, x, B1, C2);
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
	MatMultOneDim(n,A,B1,C1);					
												
	/*We perform the subtraction δx=||b-Ax'|| and place the results in the matrix C2. Then, we sort them using the qsort algorithm in order to find the maximum element in C2, which becomes the numerator of the fraction*/
	MatSubtraction(n, B2, C1, C2);
	qsort(C2, n, sizeof(double), compareDoubles);
												
	printf("\nThe absolute remainder is:\n");
    printf("\n%10.6f / %10.6f\n\n",C2[0], x[0]);
					
}












