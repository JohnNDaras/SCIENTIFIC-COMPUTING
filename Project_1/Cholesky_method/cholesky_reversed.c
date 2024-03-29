// /* A.M.:1115201800040 */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define  SIZE 25
 #define NDIM       (20)        /* number of total dimensions.  Never changes */



/*********************************************************************************************************************
*  Inversion of a symmetric matrix by Cholesky decomposition.     						     *
*  The matrix must be positive definite.                         			    			     * 
* --------------------------------------------------------------                              			     *
* 						                                                                     *
*                      												     *
*                                                        							     *
* ---------------------------------------------------------------------------------------------			     * 
* SAMPLE RUN:                                                                       				     *
*                                                                                                                    *
* Inversion of a square real symetric matrix by Cholevsky method                                                     *
* (The matrix must positive definite).                                                                               *
*
	Size of the matrix is 4:   The MATRIX is:



	|	1	 	1	 	1	 	1	 	1	 	1	 	1	 	1	 |
	|	1	 	2	 	3	 	4	 	5	 	6	 	7	 	8	 |
	|	1	 	3	 	6	 	10	 	15	 	21	 	28	 	36	 |
	|	1	 	4	 	10	 	20	 	35	 	56	 	84	 	120	 |
	|	1	 	5	 	15	 	35	 	70	 	126	 	210	 	330	 |
	|	1	 	6	 	21	 	56	 	126	 	252	 	462	 	792	 |
	|	1	 	7	 	28	 	84	 	210	 	462	 	924	 	1716 	 |
	|	1	 	8	 	36	 	120	 	330	 	792	 	1716		3432     |

Find A^(-1) with cholesky method

 Matrix A:
   1.000000   1.000000   1.000000   1.000000   1.000000   1.000000   1.000000   1.000000
   1.000000   2.000000   3.000000   4.000000   5.000000   6.000000   7.000000   8.000000
   1.000000   3.000000   6.000000  10.000000  15.000000  21.000000  28.000000  36.000000
   1.000000   4.000000  10.000000  20.000000  35.000000  56.000000  84.000000 120.000000
   1.000000   5.000000  15.000000  35.000000  70.000000 126.000000 210.000000 330.000000
   1.000000   6.000000  21.000000  56.000000 126.000000 252.000000 462.000000 792.000000
   1.000000   7.000000  28.000000  84.000000 210.000000 462.000000 924.000000 1716.000000
   1.000000   8.000000  36.000000 120.000000 330.000000 792.000000 1716.000000 3432.000000

 Matrix Inv(A):
   8.000000 -28.000000  56.000000 -70.000000  56.000000 -28.000000   8.000000  -1.000000
 -28.000000 140.000000 -322.000000 434.000000 -364.000000 188.000000 -55.000000   7.000000
  56.000000 -322.000000 812.000000 -1162.000000 1016.000000 -541.000000 162.000000 -21.000000
 -70.000000 434.000000 -1162.000000 1742.000000 -1579.000000 865.000000 -265.000000  35.000000
  56.000000 -364.000000 1016.000000 -1579.000000 1476.000000 -830.000000 260.000000 -35.000000
 -28.000000 188.000000 -541.000000 865.000000 -830.000000 478.000000 -153.000000  21.000000
   8.000000 -55.000000 162.000000 -265.000000 260.000000 -153.000000  50.000000  -7.000000
  -1.000000   7.000000 -21.000000  35.000000 -35.000000  21.000000  -7.000000   1.000000

 Here is a verification : 
 Verification A * Inv(A) = I:
   1.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
   0.000000   1.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
   0.000000   0.000000   1.000000   0.000000   0.000000   0.000000   0.000000   0.000000
   0.000000   0.000000   0.000000   1.000000   0.000000   0.000000   0.000000   0.000000
   0.000000   0.000000   0.000000   0.000000   1.000000   0.000000   0.000000   0.000000
   0.000000   0.000000   0.000000   0.000000   0.000000   1.000000   0.000000   0.000000
   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   1.000000   0.000000
   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   1.000000

To apolyto sxetiko sfalma einai:

  0.000000 /   1.000000


To apolyto sxetiko ypoloipo einai:

  0.000000 /   1.000000
***********************************************************************************************************************/



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
        for (j = i + 1; j < n; j++) {
          a[i][j] = 0;
    }
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
        for(j=0;j<n;j++)
            a[i+1][j+1]=A[i][j];

    double b[SIZE];
    for(i=0;i<n;i++)
       b[i+1]=B[i];

    for(i=1;i<=n;i++)
        a[i][n+1]=b[i];

    for(k=1;k<=n-1;k++)
    {
		if(a[k][k]==0.0)
			printf("error");
		else
		{
			pivot = a[k][k];
			for(j=k;j<=n+1;j++)
				a[k][j]= a[k][j]/pivot; 
			for(i=k+1;i<=n;i++)
			{
				factor = a[i][k];
				for(j = k;j<=n+1;j++)
					a[i][j] = a[i][j] - factor * a[k][j];
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
		// printf("%lf\n",x[i]);
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
	 
	 printf("\n\nlyseis\n");
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
*    SUBTRACTION OF TWO REAL    	  *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*1                *                                     
*            B  MATRIX N*1                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*1 = A-B          *                                     
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
*    SUBTRACTION OF TWO SQUARE REAL       *                                     
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
*                             		  *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*1        	  *                                     
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






// main program to demonstrate the use of function cholsl()
void main() {
  MAT A, A1, B, C, C4,Am1, Antistrofos={	{8.0, -28.0, 56.0, -70.0, 56.0, -28.0, 8.0, -1},
						{-28.0, 140.0, -322.0, 434.0, -364.0,  188.0, -55.0, 7.0},
						{56.0,  -322.0,   812.0,  -1162.0,  1016.0, -541.0, 162.0, -21.0},
						{-70.0,  434.0,   -1162.0,  1742.0,  -1579.0, 865.0, -265.0, 35.0},      
						{56.0, -364.0, 1016.0, -1579.0, 1476.0, -830.0, 260.0, -35.0},      
						{-28.0, 188.0, -541.0, 865.0, -830.0, 478.0, -153.0,  21.0},
						{8.0,    -55.0, 162.0, -265.0, 260.0, -153.0,  50.0,  -7.0},
						{-1.0,   7.0, 	-21.0,  35.0, -35.0,  21.0,  -7.0,   1.0} },	      
  Monadiaios={	{1.000000,   0.000000,   0.000000,   0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000,   1.000000 ,  0.000000,   0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000,   0.000000,   1.000000,   0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000,   0.000000,   0.000000 ,  1.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000,   0.000000,   0.000000 ,   0.000000, 1.000000, 0.000000, 0.000000, 0.000000},
		{0.000000,   0.000000,   0.000000 ,   0.000000, 0.000000, 1.000000, 0.000000, 0.000000},
		{0.000000,   0.000000,   0.000000 ,   0.000000, 0.000000, 0.000000, 1.000000, 0.000000},
		{0.000000,   0.000000,   0.000000 ,   0.000000, 0.000000, 0.000000, 0.000000, 1.000000}	};
  VEC p;
  int i, j, n = 8, r, flag=1, intpart, c=0, y=1, choice;
  double decpart, B1[SIZE], B2[SIZE], arr[NDIM][NDIM], C1[SIZE], C2[SIZE];
  printf(" Inversion of a square real symetric matrix by Cholesky method\n");
  printf(" (The matrix must positive def.).\n");
  printf("\nSolution of exercise 2.b.1\n");
  printf("____________________________________________\n\n\n");
  printf("                                //////////////////////////////////////////\n");
  printf("                                //      1.Read from file 2a4.txt        //\n");
  printf("                                /////////////////////////////////////////\n");               	
  printf("\n\n\nSize of the matrix is 4:   ");
  readmatrix(n, n, Am1, "2a4.txt");
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


  printf("\nFind A^(-1) with cholesky method\n");

   ////////////////////////               Find A^(-1)                      ///////////////////////////////////////////////

  char answer;
  MatPrint("Matrix A:",n,Am1);
  cholsl(n,Am1,B);
  MatPrint("Matrix Inv(A):",n,B);


  printf("\n Here is a verification : ");
  MatMult(n,A,B,C);
  MatPrint("Verification A * Inv(A) = I:",n,C);
					  

					  		                            /*Calculation of Errors.*/		
					
// i]						
   /*We perform the subtraction δx=||A^(-1) - A'^(-1)|| and place the results in matrix C4. Then, we find the sums of the columns for each row and store them in matrix C2. */
  MatSubtraction2(n, B, Antistrofos, C4);
  MatRowsSum(n, C4, C2);
				  
   /*Afterwards, I sort them using the qsort algorithm in order to find the maximum element in C2, which becomes the numerator of the fraction*/
  qsort(C2, n, sizeof(double), compareDoubles);									
				
   /* We find the sums of the columns for each row of the inverse we have found using the Cholesky method and store them in matrix C1. Then, we sort the elements of C1 using the qsort algorithm in order to find the maximum element, which becomes the denominator of the fraction.*/
   MatRowsSum(n, B, C1);
   qsort(C1, n, sizeof(double), compareDoubles);
	  

   printf("\nThe absolute relative error is: \n");		
   printf("\n%10.6f / %10.6f\n\n",C2[0], C1[0]);


//ii]						
	/*I perform the subtraction δx=||A*A^(-1) - I || and place the results in matrix C4. Then, I find the sums of the columns for each row and store them in matrix C2.*/
   MatSubtraction2(n, C, Monadiaios, C4);
   MatRowsSum(n, C4, C2);			  
				 
   /* Afterward, I sort them using the qsort algorithm in order to find the maximum element in C2, which becomes the numerator of the fraction*/					 
   qsort(C2, n, sizeof(double), compareDoubles);

   /* We find the sums of the columns for each row of the inverse that I have found using the Cholesky method and store them in matrix C1. Then, we sort the elements of C1 using the qsort algorithm to find the maximum element, which becomes the denominator of the fraction */				
   MatRowsSum(n, B, C1);
   qsort(C1, n, sizeof(double), compareDoubles);
   printf("\nThe absolute relative remainder is:\n");		
   printf("\n%10.6f / %10.6f\n\n",C2[0], C1[0]);					  
					  
}
