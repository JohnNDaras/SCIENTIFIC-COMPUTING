// ESOR
// 1115201800040
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

double ** Matrix_5x5(int *, double**);
double ** Insert_A(int*);
double * Insert_X(int);
int GivenMatrix();
double * Insert_rand_X(int);
double ** Read_file_Axb(int*, double**);
void ESOR_mainmethod (double**, double*, int);
double ** dispArray(int, int);
double sign();
void ESOR_iterativeMethod (double **, double **, int, double, double, int, int *);
double ** Matrix_10x10(int*, double**);
double ** Matrix_NxN (int*, double**);
double ** Insert_rand_A(int*);

struct timeval t1, t2, t3, t4;

int main(int argc, char* argv[]) {
	int back2menu, n,Selection;
	double **A, *X;
	printf("\nLINEAR SYSTEM RESOLUTION OF Ax = b SYSTEMS USING ESOR ITERATIVE METHOD\n");
	while(1) {
	back2menu = 0;

	//main menu
	do{
	printf ("\n\nMAIN MENU\n\n\t1. Please enter the pentadiagonal matrix A and the vector b from the keyboard.\n\t2. Choose from a list of provided matrices.\n\t3. Select the dimensions of a random matrix.\n\t4. Enter a matrix from a file named 'ask2_ESOR.txt'.\n\n\t5. Exit.\n");
			printf("INPUT:  ");
			scanf("%d", &Selection);
	} while((Selection<1)||(Selection>5));

	switch(Selection) {
	case 1: A = Insert_A(&n);
		X = Insert_X(n);
		break;
	case 2: switch(GivenMatrix()) {
				case 1: A = Matrix_5x5(&n, &X);
					break;
				case 2: A =Matrix_10x10(&n, &X);
					break;
				case 3: A = Matrix_NxN (&n, &X);
					break;
				case 4: back2menu = 1;
					break;
				default: break;
						}
				break;
	case 3: A = Insert_rand_A(&n);
				X = Insert_rand_X(n);
				break;
	case 4: A = Read_file_Axb(&n, &X);
				break;


	case 5: 
		return 0;
	default: break;
	}
	if (back2menu!=1)
		ESOR_mainmethod(A, X, n);
	}
}



void ESOR_iterativeMethod (double **A, double **x1, int n, double t, double w, int itmax, int *itcount) {
	//algorithm of ESOR for factorizing matrix A in L and U using partial pivoting
	long double e = 0.0000005;
	double *x0, *k, *l, *d, *r, *s, approximation;
	int i, j;
	//xo, k, l, d, r, s matrix allocations
	x0 = (double*)malloc(n*sizeof(double));
	k = (double*)malloc(n*sizeof(double));
	l = (double*) malloc(n*sizeof(double));
	d = (double*)malloc(n*sizeof(double));
	r = (double*)malloc(n*sizeof(double));
	s = (double*) malloc(n*sizeof(double));
	//initialization of xo, k, l, r, S
	for(i=0; i<n; i++) {
	x0[i] = 0;
	k[i] = 0;
	l[i] = 0;
	d[i] = 0;
	r[i] = 0;
	s[i] = 0;
	}
	//x0, k, l, d, r, s assignements
	for(i=0; i<n; i++){
	    x0[i] = A[i][n];
	    for(j=0; j<n; j++)
		switch(abs(i-j)){
				case 0: d[i] = A[i][j];
					break;
				case 1: 
					if(i>j) l[i]= A[i][j];
					else r[i]=A[i][j];
					break;
				case 2:
					if(i>j) k[i]=A[i][j];
					else s[i]=A[i][j];
					break;
				default:	
					break;
		}
	}


	while((*itcount) < itmax) {
		approximation = 0.0;
		//1 - compute values of x1
		for(i=0; i<n; i++)
		(*x1)[i] = (1.0-t)*x0[i] - w* ((k[i]/d[i])*(*x1)[i-2]+(l[i]/d[i])*(*x1)[i-1]) -
						(t-w)*((k[i]/d[i])*x0[i-2]+(l[i]/d[i])*x0[i-1]) -
						t*((r[i]/d[i]) *x0[i+1]+(s[i]/d[i])*x0[i+2]-(A[i][n]/d[i]));
	
		//2 increment of iterations' count
		(*itcount)++;
	
		for(i=0; i<n; i++)
		if(fabs((*x1)[i]-x0[i])>approximation)
				approximation=fabs((*x1)[i]-x0[i]);
		if(approximation <e)
				//3 - check if approximation is less than desired accuracy and if yes stop method
				break;
		else
		//4 - assign x1 to xo to aprroximate better in next iteration
		for(i=0; i<n; i++)
		x0[i] = (*x1)[i];
	}
	free (x0);
	free(k);
	free(l);
	free(d);
	free(r);
	free(s);
}


void ESOR_mainmethod (double **A, double *x, int n) {

	double  t_opt, w_opt,t, w, cpuduration, itcpuduration, optimal_cputime;
	int i, j, itmax, min_itcount=0,itcount=0,Select;



	printf("Enter the maximum number of iterations you wish to be applied for computing the solution of the system:\n  ");
	scanf("%d", &itmax);


	do{
		printf ("IMPLEMENTATION MENU\n\t1. Experimental validation of the accuracy of the ESOR algorithm for the analysis of the linear system Ax=b with x=(1,1, ...,1)^T\n\t2. Experimental study of the convergence of the ESOR iterative method with parameters t,w belonging to [0.1, 1.9] with a step of 0.1\n\t3. Exit. \n");
		printf("INPUT: ");
		scanf("%d", &Select);
	}while( (Select<1)||(Select>3));

	switch(Select) {
		case 1: printf("Enter the parameters t and w, where t, w belong to [0.1, 1.9] (e.g., t=0.1, w=0.4):\n");
			printf("t = ");
			scanf("%lf", &t);
			printf("w = ");
			scanf("%lf", &w);
			gettimeofday(&t1, NULL);
			ESOR_iterativeMethod(A ,&x ,n ,t,w ,itmax ,&itcount);
			gettimeofday(&t2, NULL);
			cpuduration = (t2. tv_sec - t1.tv_sec) * 1000.0;
			cpuduration += (t2. tv_usec - t1.tv_usec) / 1000.0;

			(itcount < itmax) ? printf("\nThe approximate solution for the pentadiagonal matrix was found with parameters t=%2.1lf and w=%2.1lf after %d iterations.n", t,w,itcount) :
								printf("\nThe approximate solution for the pentadiagonal matrix was not found within %d iterations with parameters t=%2.1lf and w=%2.1lf. The latest approximate solution is as follows.\n",itmax, t, w);


			for(i=0; i<n; i++, printf("\n"))
				for (j=0; j<n+1 ; j++)
					printf("%7.1lf", A[i][j]);
			printf("\n\nThe solution of the system is: \n");
			printf("\tx = ( ");
			for(i=0; i<n-1; i++)
				printf("%7.3lf, ", x[i]);
			printf("%7.3lf )\n\n", x[n-1]);

			printf("\nThe execution time is: %.41f ms\n\n", cpuduration);
			break;

		 case 2:  
			printf("The parameters t and w are iteratively selected within the range [0.1, 1.9] in order to find the optimal values for t and w.\n");
			printf("__________________________________________________________________________________________________________________________\n");
			printf("\n////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
			for(t=0.1 ; t<2.0 ; t+=0.1)
				for(w=0.1 ; w<2.0 ; w+=0.1){
					itcpuduration = 0.0;
					gettimeofday(&t1, NULL);
					ESOR_iterativeMethod(A, &x, n, t, w, itmax, &itcount);
					gettimeofday(&t2, NULL);
					itcpuduration = (t2. tv_sec - t1.tv_sec) * 1000.0;
					itcpuduration += (t2.tv_usec - t1.tv_usec) / 1000.0;
					cpuduration += itcpuduration;
					//printing each time only when a better solution is found
					if((itcount < min_itcount) || (t==0.1 && w==0.1)) {
					min_itcount = itcount;
					t_opt = t;
					w_opt = w;
					printf("//   The best parameters up to these iterations are t_optimal = %2.1lf and w_optimal = %2.1lf  //\n", t_opt, w_opt);
					}
					itcount = 0;
				}
			printf("////////////////////////////////////////////////////////////////////////////////////////////////////////\n");

			if (min_itcount < itmax) 
				printf("\nThe approximate value of the solution for the matrix was found with optimal parameters t=%2.11F and w=%2.11f after %d iterations. \n", t_opt, w_opt, min_itcount) ;
			else    
				printf("\nThe approximate value of the solution for the matrix was not found within %d iterations with parameters t=%2.1lf and w=%2.1lf. The latest approximate solution is as follows.\n", itmax, t_opt, w_opt);

			gettimeofday(&t3, NULL);
			ESOR_iterativeMethod(A, &x, n, t_opt, w_opt, itmax, &min_itcount);
			gettimeofday(&t4, NULL);
			optimal_cputime = (t4.tv_sec - t3.tv_sec) * 1000.0;
			optimal_cputime += (t4.tv_usec - t3.tv_usec) / 1000.0;
			
			for(i=0; i<n; i++, printf("\n"))
				for (j=0; j<n+1 ; j++)
					printf("%7.11f", A[i][j]);

			printf("\n\nThe solution of the system is: \n");
			printf("\tx = ( ");
			for(i=0; i<n-1; i++)
			printf("%7.3lf, ", x[i]);
			printf("%7.3lf)\n\n", x[n-1]);
			printf("The optimal parameter t is: %2.1lf\n", t_opt);
			printf("The optimal parameter w is:%2.1lf\n", w_opt);
			printf("The number of iterations executed with the optimal parameters t and w is: %d\n", min_itcount);
			printf("\nThe execution of the PSD method for the linear system Ax=b with all values of t and w within the range took: %.41f ms\n\n", cpuduration);
			printf("\nThe execution of the PSD with the minimum number of iterations took:%.41f ms\n\n", optimal_cputime);
			break;
	exit(0);
	}

	for(i=0; i<n; i++)
		free (A[i]);
	free(A);
	free(x);
}



double ** Insert_A(int *n) {
	//data of matrix A and vector b are inserted by the user
	double **A;
	int i, j;
	printf("Please enter the dimension N of the matrix: ");
	scanf("%d", n);
	A = dispArray(*n, *n+1); //the one more column represents the b vector. this allowance has effect in the whole written program
	printf("\nEnter below the elements of the matrix:\n");
	for(i=0; i<*n; i++) {
		for(j=0; j<*n; j++){
			printf("A[%d] [%d] = ", i, j);
			scanf ("%lf", &A[i][j]);
		}
	}

	printf("\nEnter the data for b: \n");
	for (i=0; i<*n; i++){
		printf("b[%d] = ", i);
		scanf ("%lf", &A[i][*n]);
	}
	return A;
}



double * Insert_X(int n) {
	//data of matrix X are inserted by the user
	double *X;
	int i;
	X = (double*)malloc(n*sizeof(double));
	printf("\nPlease enter the vector X (e.g., (1,1,1,1,1)^T for a 5x5 pentadiagonal matrix): \n");
	for (i=0; i<n; i++) {
		printf("X[%d] = ", i);
		scanf("%lf", &X[i]);
	}
	return X;
}



int GivenMatrix() {
	int matrix_no;
	do {
	printf("\n\nSpecific matrices LIST		\nSelect a specific matrix (from the given data in the statement).:\n\t1.Pentadiagonal matrix with dimensions 5X5 (adjustment 1.1).\n\t2. Pentadiagonal matrix with dimensions 10X10 (adjustment 1.ii)\n\t3. Pentadiagonal matrix with dimensions NXN, where N, a, b, c, and d are given by you with values 100, 1000, 10000 (for N) and 0.1-1.9 (for each of a, b, c, d)\n\n\t4. Back to menu.\n");
	printf("INPUT : ");
	scanf("%d", &matrix_no);
	}while((matrix_no<1)||(matrix_no>4));
	return matrix_no;
}



double ** Matrix_5x5(int *n, double **X) {
	//returns a specific matrix (adjustment 1.i)
	double **A;
	*n=5;
	A = dispArray(*n, *n+1);          //the one more column corresponds to the b vector
	A[0][0]=4; 			A[0][1]=-0.3; 		A[0][2]=-0.4; 		 A[0][3]=0;	 	A[0][4]=0; 		A[0][5]=3.3;
	A[1][0]=-0.2; 			A[1][1]=4; 		A[1][2]=-0.3; 		 A[1][3]=-0.4;		A[1][4]=0; 		A[1][5]=3.1;
	A[2][0]=-0.1; 			A[2][1]=-0.2; 		A[2][2]=4;		 A[2][3]=-0.3; 		A[2][4]=-0.4; 		A[2][5]=3;
	A[3][0]=0; 			A[3][1]=-0.1; 		A[3][2]=-0.2; 		 A[3][3 ]=4; 		A[3][4]=-0.3; 		A[3][5]=3.4;
	A[4][0]=0;			A[4][1]=0; 		A[4][2]=-0.1;		 A[4][3]=-0.2; 		A[4][4]=4; 		A[4][5]=3.7;

	*X = (double*)malloc(*n*sizeof(double));
	(*X)[0]=1; 		(*X)[1]=1; 		(*X) [2]=1;		 (*X) [3]=1; 		(*X) [4]=1;
	return A;
}


double ** Matrix_10x10(int *n, double **X) {
	//returns a specific matrix (adjustment 1.ii)
	double **A;
	*n=10;
	A = dispArray(*n, *n+1);		//the one more column corresponds to the b vector

	A[0][0]=4; 			A[0][1]=-0.2; 		A[0][2]=-0.1; 		A[0] [3]=0; 		A[0][4]=0; 		A[0][5]=0; 		A[0][6]=0; 		A[0][7]=0; 		A[0][8]=0; 		A[0][9]=0;		A[0][10]=3.7;
	A[1][0]=-0.3;			A[1][1]=4; 		A[1][2]=-0.2; 		A[1][3]=-0.1; 		A[1][4]=0;		A[1][5]=0;		A[1][6]=0; 		A[1][7]=0;		A[1][8]=0;		A[1][9]=0; 		A[1][10]=3.4;
	A[2][0]=-0.4;		    	A[2][1]=-0.3; 		A[2][2]= 4;		A[2][3]=-0.2;		A[2][4]=-0.1;		A[2][5]=0;		A[2][6]=0;		A[2][7]=0;		A[2][8]=0; 		A[2][9]=0; 		A[2][10]=3;						
	A[3][0]=0;			A[3][1]=-0.4; 		A[3][2]=-0.3;		A[3][3]=4; 		A[3][4]=-0.2; 		A[3][5]=-0.1;		A[3][6]=0;		A[3][7]=0; 		A[3][8]=0;		A[3][9]=0;		A[3][10]=3;
	A[4][0]=0;			A[4][1]=0; 		A[4][2]=-0.4;		A[4][3]=-0.3; 		A[4][4]=4;		A[4][5]=-0.2;		A[4][6]=-0.1; 		A[4][7]=0; 		A[4][8]=0;		A[4][9]=0;		A[4][10]=3;
	A[5][0]=0; 			A[5][1]=0; 		A[5][2]=0; 		A[5][3]=-0,4;		A[5][4]=-0,3;		A[5][5]=4;		A[5][6]=-0,2;		A[5][7]=-0.1; 		A[5][8]=0; 		A[5][9]=0;		A[5][10]=3;					
	A[6][0]=0;			A[6][1]=0;		A[6] [2]=0;		A[6][3]=0;		A[6][4]=-0.4; 		A[6][5]=-0.3;		A[6][6]=4; 		A[6][7]=-0.2; 		A[6][8]=-0.1; 		A[6][9]=0;		A[6][10]=3;
	A[7][0]=0;			A[7][1]=0;		 A[7][2]=0;		A[7][3]=0;		A[7][4]=0; 		A[7][5]=-0.4; 		A[7][6]=-0.3; 		A[7][7]=4; 		A[7][8]=-0.2; 		A[7][9]=-0.1; 		A[7][10]=3;
	A[8][0]=0;			A[8][1]=0;		A[8][2]=0; 		A[8][3]=0;		A[8][4]=0; 		A[8][5]=0; 		A[8][6]=-0.4; 		A[8][7]=-0.3; 		A[8][8]=4; 		A[8][9]=-0.2;		A[8][10]=3.1;
	A[9][0]=0;			A[9][1]=0; 		A[9] [2]=0; 		A[9][3]=0; 		A[9][4]=0; 		A[9][5]=0; 		A[9][6]=0; 		A[9][7]=-0.4; 		A[9][8]=-0.3; 		A[9][9]=4;		A[9][10]=3.3;

	*X = (double*) malloc(*n*sizeof(double));
	(*X)[0]=1;		 (*X) [1]=1; 		(*X) [2]=1; 		(*X) [3]=1; 		(*X) [4]=1; 		(*X) [5]=1; 		(*X) [6]=1;		 (*X) [7]=1; 		(*X) [8]=1; 		(*X) [9]=1;

	return A;
}




double ** Matrix_NxN (int *n, double **X) {
	//returns a specific matrix (adjustment 2)
	double **A, a, b, c, d;				//where -a=k, -b=l, -c=r, -d=s
	int i, j;
	printf("Please enter the dimension N of the pentadiagonal matrix.[NXN]  (p.x. 100, 1000, 10000) : ");
	scanf("%d", n);
	printf("Please enter the parameters a, b, c, d in that order, with values for a, b, c, d in the range 0.1(0.1)1.9.(i.e. a=1.2, b=0.9, c=0.6, d=0.3) :\n");
	printf("a = ");
	scanf("%lf", &a);
	printf("b=  ");
	scanf("%lf",&b);
	printf("c= ");
	scanf("%lf", &c);
	printf("d = ");
	scanf("%lf", &d);

	*X = (double*) malloc(*n*sizeof(double));
	A = dispArray(*n, *n+1);		//the one more column corresponds to the b vector
	//values for matrix A
	for(i=0; i<*n; i++)
		for(j=0; j<*n; j++)
			switch(abs(i-j)){
				case 0: A[i][j]=4;
					break;
				case 1: 
					if(i>j) A[i][j]=-b;
					else A[i][j]=-c;
					break;
				case 2: 
					if(i>j) A[i][j]=-a;
					else A[i][j]=-d;
					break;
				default:A[i][j]=0.0;
					break;

	}

	//values for X
	for(i=0; i<*n; i++)
	(*X)[i]=1;
	//b computation
	for(i=0; i<*n; i++) {
		A[i][*n] = 0;		//initializes into zero all the values of b vector
		for(j=0; j<*n; j++)
			A[i][*n] += A[i][j]*(*X) [j];
		//the above vector b cause of x=(1,1, ...,1)^T means b[i] = a[i]+b[i]+diag[i]+c[i]+d[i] = -(k[i]+2[i]+d[i]+r[i]+s[i])
	}
	return A;
}



double ** Insert_rand_A(int *n) {
	srand(time(NULL));
	double **A, k, l, r, s, values, diag_value, mult; 
	int i, j, case_values;
	printf("Please enter the dimension N of the pentadiagonal matrix [NXN] (p.x. 100, 1000, 10000) : ");
	scanf("%d", n);
	
	
	
	//values for parameters of matrix A
	values = 0.1+(0.1*(double) (rand()%19));	//desirable values of parameters a,b,c,d between [0.1, 1.9]
	case_values = (int) (values*10);
	diag_value = (sign()*(4.0+(double)((rand()%385)/4)))-1;		//desirable values of diagonial between [-100,-4] or [4,100] all multiplies of 4(-100, ...,-4,4,8,12, ...,
	
	
	switch(case_values) {
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
			mult = 0.1*(1.0+(double)(rand()%4));
			k = sign() *values;
			l = sign()*(k+mult);
			r = sign()*(k+2.0*mult);
			s = sign()*(k+3.0*mult);
			break;
		case 6:
		case 7:
		case 8:
		case 9:
		case 10:
			mult = 0.1*(1.0+(double)(rand()%3));
			k = sign() *values;
			l = sign()*(k+mult);
			r = sign()*(k+2.0*mult);
			s = sign()*(k+3.0*mult);
			break;
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
			mult = 0.1*(1.0+(double)(rand()%3));
			s = sign() *values;
			r = sign()*(s-mult);
			l = sign()*(s-2.0*mult);
			k = sign()*(s-3.0*mult);
			break;
		case 16:
		case 17:
		case 18:
		case 19:
			mult = 0.1*(1.0+(double) (rand()%4));
			s = sign() *values;
			r = sign()*(s-mult);
			l = sign()*(s-2.0*mult);
			k = sign()*(s-3.0*mult);
			break;
		default:break;
	}
	printf("\n\n\nPentadiagonial matrix parameters selected as follows:\nk=%2.11f, L=%2.11f, d=%2.11f, r=%2.11f, s=%2.11f\n\n\n", k, l, diag_value,r,s);
	A = dispArray(*n, *n+1);			//the one more column corresponds to the b vector



	//values for matrix A
	for(i=0; i<*n; i++)
		for(j=0; j<*n; j++)
			switch(abs(i-j)){
			case 0: A[i][j]=diag_value;
				break;
			case 1: 
				if(i>j) A[i][j]=l;
				else A[i][j]=r;
				break;
			case 2: 
				if(i>j) A[i][j]=k;
				else A[i][j]=s;
				break;
			default:A[i][j]=0.0;
				break;
		}

	//randomly selected data for b vector below so that X is (1,1, ...,1,1)^T
	for(i=0; i<*n; i++) {
		A[i][*n]=0.0;
		//initializes into zero all the values of b vector
		for(j=0; j<*n; j++)
		A[i][*n] += A[i][j];
		//the above vector b cause of x=(1,1, ...,1)^T means b[i] = a[i]+b[i]+diag[i]+[i]+d[i] = -(k[i]+2[i]+d[i]+r[i]+s[i])
	}
	return A;
}




double * Insert_rand_X(int n) {
	//data of matrix X are randomly selected and inserted (adjustment 3)
	double *X;
	int i;
	X = (double*) malloc(n*sizeof(double));
	for (i=0; i<n; i++)
	X[i] = 1;
	return X;
}



double ** Read_file_Axb(int *n, double **X) {
	double **A;
	int i, j;
	FILE *file;
	if((file=fopen("ask2_ESOR.txt", "r"))==NULL) {
			perror("FILE ERROR : Either the file with name 'ask2_ESOR.txt' does not exist or it is placed in a wrong directory");
			exit(0);
	}

	fscanf(file, "%d", n);				
	printf("Dimensions given from the file are %d/%d\n", *n, *n);

	A = dispArray(*n, *n+1);
	//the one more column corresponds to the b vector
	(*X) = (double*)malloc(*n*sizeof(double));

	//reading data from file and import them in matrix A below
	for(i=0; i<*n; i++)
		for(j=0; j<*n; j++)
			fscanf(file, "%lf",&A[i][i]);

	//reading data from file and import them in matrix b below
	for(i=0; i<*n; i++)
		fscanf(file, "%lf",&A[i][*n]);

	//reading data from file and import them in vector X below
	for(i=0; i<*n; i++)
		fscanf(file, "%lf",&(*X)[i]);

	fclose(file);
	return A;
}



double ** dispArray(int n, int m) {
//dynamic allocation of matrices
	double **A;
	int i, j;
	A= (double**)malloc(n*sizeof(double*));
	for(i=0;i<n;i++)
		*(A+i)=(double*) malloc(m*sizeof(double));
	
	return A;
}



double sign() {
	if (rand()%2==0)
	return 1.0;
	else
	return -1.0;
}


