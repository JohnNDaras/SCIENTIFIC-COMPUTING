// PSD
// 1115201800040

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

int menu();
double **Insert_A(int*);
double * Insert_X(int);
int GivenMatrix();
double **Matrix_5x5(int *, double**);
double **Matrix_10x10(int*, double**);
double ** Matrix_NxN (int*, double**);
double ** Read_file_Axb(int*, double**);
int PSD_menu();
void PSD_mainmethod (double**, double*, int);
double ** dispArray(int, int);
double sign();
double ** Insert_rand_A(int*);
double * Insert_rand_X(int);
void PSD_iterativeMethod (double **, double **, int, double, double, int, int *);
struct timeval t1, t2, t3, t4;




int main(int argc, char* argv[]) {
    int b, n, Select;
    double **A, *X;

    while(1) {
        b = 0;

        do
	{
            printf ("\n\nMAIN MENU\n\n\t1. Parakalw eisagete ton pentadiagwnio pinaka A kai to dianysma b apo to plhkrologio\n\t2. Epilekste apo mia lista apo dosmenoys pinakes\n\t3. Epilekste thn diastash enos tyxaioy pinaka\n\t4. Eisagwgh pinaka apoarxeio me onoma: 'ask2_ESOR.txt'\n\n\t5. Eksodos.\n");	
            printf("INPUT: ");
            scanf("%d", &Select);
        }while( (Select<1)||(Select>5));

        switch(Select) {
        case 1: A = Insert_A(&n);
		X = Insert_X(n);
		break;
        case 2: switch(GivenMatrix()) {
			case 1: A = Matrix_5x5(&n, &X);
						break;
			case 2: A = Matrix_10x10(&n, &X);
						break;
			case 3: A = Matrix_NxN(&n, &X);
						break;
			case 4: b = 1;
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

    	default:break;
    }
    if (b!=1)
    	PSD_mainmethod(A, X, n);
    }
}



void PSD_iterativeMethod (double **A, double **x, int n, double t, double w, int itmax, int *itcount) {
    //algorithm of PSD for factorizing matrix A in L and U using partial pivoting
    long double e = 0.0000005;
    double *x0, *x1, *Ux0, *Ux1, *k, *l, *d, *r, *s, approximation;
    int i, j;


    //xo, x1, Uxo, Ux1, k, l, d, r, s matrix allocations
    x0 = (double*)malloc(n*sizeof(double));
    x1 = (double*)malloc(n*sizeof(double));
    Ux0 = (double*)malloc(n*sizeof(double));
    Ux1 = (double*)malloc(n*sizeof(double));
    k = (double*) malloc(n*sizeof(double));
    l = (double*)malloc(n*sizeof(double));
    d = (double*)malloc(n*sizeof(double));
    r = (double*)malloc(n*sizeof(double));
    s = (double*)malloc(n*sizeof(double));
    //initialization of xo, x1, Uxo, Ux1, k, l, r, s
    for(i=0; i<n; i++) {
        x0[i] = 0;
        x1[i] = 0;
	Ux0[i] = 0;
	Ux1[i] = 0;
	k[i]= 0;
	l[i] = 0;
	d[i] = 0;
	r[i] = 0;
	s[i] = 0;
    }

    for(i=0; i<n; i++){
	x0[i] = A[i][n];
	x1[i] = 0.0;
	Ux0[i] = 0.0;
	Ux1[i] = 0.0;
	for(j=0; j<n; j++)
	    switch(abs(i-j)){
		case 0: d[i] = A[i][j];
			break;
		case 1: if (i>j) 
			    l[i]=A[i][j] ;
			else 
			    r[i]=A[i][j];
			break;
		case 2: if (i>j) 
			    k[i]=A[i][j] ; 
			else 
			    s[i]=A[i][j];
			break;
		default: break;
		}
    }
    while((*itcount) < itmax) {
	approximation = 0.0;
	//1 compute values of x
	for(i=0; i<n; i++){
		if((*itcount) == 0)
			Ux0[i] = -((r[i]/d[i])*x0[i+1]+(s[i]/d[i])*x0[i+2]);
		x1[i] = (1.0-t)*x0[i] - w*((k[i]/d[i])*x1[i-2]+(l[i]/d[i])*x1[i-1]) -
			(t-w)*((k[i]/d[i])*x0[i-2]+(l[i] / d[i])*x0[i-1]) -
		t*((r[i]/d[i]) *x0[i+1]+(s[i]/d[i])*x0[i+2]-(A[i][n]/d[i]));
		Ux1[i] = -((r[i]/d[i])*(*x)[i+1]+(s[i]/d[i])*(*x)[i+2]);
		(*x)[i] = x1[i]+w* (Ux1[i] -Ux0[i]);
		Ux0[i] = Ux1[i];
	}
	//increment of iterations' count
	//2
	(*itcount)++;
	for(i=0; i<n; i++)
	approximation = 0.0;
	//1 compute values of x
	for(i=0; i<n; i++){
	    if((*itcount) == 0)
		Ux0[i] = -((r[i]/d[i])*x0[i+1]+(s[i]/d[i])*x0[i+2]);
	    x1[i] = (1.0-t)*x0[i] - w*((k[i]/d[i])*x1[i-2]+(l[i]/d[i])*x1[i-1]) -
			(t-w)*((k[i]/d[i])*x0[i-2]+(l[i] / d[i])*x0[i-1]) -
	    t*((r[i]/d[i]) *x0[i+1]+(s[i]/d[i])*x0[i+2]-(A[i][n]/d[i]));
	    Ux1[i] = -((r[i]/d[i])*(*x)[i+1]+(s[i]/d[i])*(*x)[i+2]);
	    (*x)[i] = x1[i]+w* (Ux1[i] -Ux0[i]);
	    Ux0[i] = Ux1[i];
	}
	//increment of iterations' count
	//2
	(*itcount)++;
	for(i=0; i<n; i++)
	    if(fabs((*x)[i]-x0[i])>approximation)
	        approximation=fabs((*x)[i]-x0[i]);
	if(approximation <e)
	    //3 check if approximation is less than desired accuracy and if yes stop method
	    break;
	else
	    //4 - assign x to xo to aprroximate better in next iteration
	    for(i=0; i<n; i++)
	        x0[i] = (*x)[i];if(fabs((*x)[i]-x0[i])>approximation)
	approximation=fabs((*x)[i]-x0[i]);
	if(approximation <e)
	    //3 check if approximation is less than desired accuracy and if yes stop method
	    break;
	else
	    //4 - assign x to xo to aprroximate better in next iteration
	    for(i=0; i<n; i++)
	        x0[i] = (*x)[i];

    }

    free (x0);
    free (x1);
    free (Ux0);
    free (Ux1);
    free(k);
    free(l);
    free(d);
    free(r);
    free(s);
}


void PSD_mainmethod (double **A, double *x, int n) {

    double t_opt, w_opt,t, w, cpuduration, itcpuduration, optimal_cputime;
    int Select,min_itcount=0,i, j, itmax,itcount=0;

    printf("Dwste ton megisto arithmo epanalhpsewn poy epiyhymeite na efarmostoyne gia ton ypologismo ths lyshs toy systhmatos:\n  ");

    scanf("%d", &itmax);

    do{
        printf ("IMPLEMENTATION MENU\n\t1. Peiramatikh epalhtheysh ths orthothtas toy PSD algorithmoy gia thn analysh toy grammikoy systhmatos Ax=b me x=(1,1, ...,1)^T\n\t2. Peiramatikh meleth ths sygklishs ths PSD epanalhptikhs methodoy me t,w na anhkoune sto [0.1, 1.9] me vhma 0.1\n\t3. Eksodos. \n");
        printf("Epilogh: ");
        scanf("%d", &Select);
    }while( (Select>3)||(Select<1));


    switch(Select) {
        case 1: printf("Dwste tis t kai w parametroys s' aythn thn seira opoy t, w anhkoyn 0.100.1)1.9 (p.x. t=0.1, w=0.4) :\n");
                printf("t = ");
		scanf("%lf", &t);
		printf("w = ");
		scanf("%lf", &w);
		gettimeofday(&t1, NULL);
		PSD_iterativeMethod(A, &x, n, t, w, itmax, &itcount);
		gettimeofday(&t2, NULL);
		cpuduration = (t2.tv_sec - t1.tv_sec) * 1000.0;
		cpuduration += (t2. tv_usec - t1.tv_usec) / 1000.0;

		if (itcount < itmax) 
		    printf("\nH proseggistikh timh ths lyshs gia ton pinaka vrethike me parametroys tis t=%2.11F kai w=%2.11f meta apo %d epanalipseis \n", t,w, itcount) ;
		else 
		    printf("\nH proseggistiki timh ths lyshs gia ton pinaka den vrethike se %d epanalhpseis me parametroys tis t=%2.1lf kai w=%2.1lf. H teleytaia proseggistikh lysh einai h parakatw. \n", itmax , t, w); 


		for(i=0; i<n; i++)
		    printf("\n");
		    for (i=0; j<n+1 ; j++)
		        printf("%7.1lf", A[i][j]);
		printf("\n\nH lysh toy systhmatos einai : \n");
		printf("\tx = ( ");
		for(i=0; i<n-1 ; i++)
		    printf("%6.4lf, ", x[i]);
		printf("%6.4lf )\n\n", x[n-1]);
		printf("\nO xronos ekteleshs einai : %.5lf ms\n\n", cpuduration);
		break;
	case 2: 
		printf("\n\n\nOi parametroi t kai w epilegontai epanalhptika sto diasthma [0.1, 1.9] wste na vroyme tis veltistes times toy t kai toy w\n");
		printf("__________________________________________________________________________________________________________________________\n");
		printf("\n////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
		
		for(t=0.1 ; t<2.0 ; t+=0.1)
		    for(w=0.1 ; w<2.0 ; w+=0.1) {
			itcpuduration = 0.0;
			gettimeofday(&t1, NULL);
			PSD_iterativeMethod(A, &x, n, t, w, itmax, &itcount);
			gettimeofday(&t2, NULL);
			itcpuduration = (t2.tv_sec - t1.tv_sec) * 1000.0;
			itcpuduration += (t2. tv_usec - t1.tv_usec) / 1000.0;
			cpuduration += itcpuduration;
			//printing each time only when a better solution is found
			if((itcount < min_itcount) || (t==0.1 && w==0.1)) {
			    min_itcount = itcount;
			    t_opt = t;
			    w_opt = w;
			    printf("//   Oi kalyteres parametroi mexri aytes tis epanalhpseis einai t_optimal = %2.1lf and w_optimal = %2.1lf   //\n", t_opt, w_opt);
			}
			itcount = 0;
		    }
		printf("////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
		if (min_itcount < itmax) 
		    printf("\nH proseggistikh timh ths lyshs gia ton pinaka vrethike me veltistes parametroys tis t=%2.11F kai w=%2.11f meta apo %d epanalipseis \n", t_opt, w_opt, min_itcount) ;
		else    
		    printf("\nH proseggistiki timh ths lyshs gia ton pinaka den vrethike se %d epanalhpseis me parametroys tis t=%2.1lf kai w=%2.1lf. H teleytaia proseggistikh lysh einai h parakatw.\n", itmax, t_opt, w_opt);
		
		gettimeofday(&t3, NULL);
		PSD_iterativeMethod(A, &x, n, t_opt, w_opt, itmax, &min_itcount);
		gettimeofday(&t4, NULL);
		optimal_cputime = (t4.tv_sec - t3.tv_sec) * 1000.0;
		optimal_cputime += (t4.tv_usec - t3.tv_usec) / 1000.0;
		
		printf("\nMatrix A with vector b in last column : \n");
		for(i=0; i<n; i++, printf("\n"))
		    for (j=0; j<n+1 ; j++)
			printf("%7.1lf", A[i][j]);
		
		printf("\n\nH lysh toy systhmatos einai : \n");
		printf("\tx = ( ");
		for(i=0; i<n-1; i++)
		    printf("%7.31f, ", x[i]);
		printf("%7.3lf )\n\n", x[n-1]);
		printf("H veltisth parametros t einai : %2.1lf\n", t_opt);
		printf("H veltisth parametros w einai : %2.1lf\n", w_opt);
		printf("O arithmos twn epanalhpsewn poy ektelesthkan me tis veltistes parametroys t kai w einai : %d\n", min_itcount);
		printf("\nH ektelesh ths PSD methodoy gia to grammiko systhma Ax=b me oles tis times t kai w sto diasthma phre : %.41f ms\n\n", cpuduration);
		printf("\nH ektelesh ths PSD me ton elaxisto arithmo epanalhpsewn phre : %.41f ms\n\n", optimal_cputime);
		break;


	case 3:  
	    exit(0);
	//not return 0..
    }
    for(i=0; i<n; i++)
        free (A[i]);
    free(A);
    free(x);
}




double **Insert_A(int *n) {
    //data of matrix A and vector b are inserted by the user
    double **Matrix;
    int i, j;
    printf("Parakalw eisagete thn diastash N toy pinaka : ");
    scanf("%d", n);
    Matrix = dispArray(*n, *n+1);
    //the one more column represents the b vector. this allowance has effect in the whole written program
    printf("\nEisagete parakatw ta stoixeia toy pinaka: \n");
    for(i=0; i<*n; i++) {
	for(j=0; j<*n; j++){
	printf("A[%d][%d] = ", i, j);
	scanf ("%lf", &Matrix[i][j]);
	}
    }
    printf("\nEisagete ta dedomena toy b: \n");
    for (i=0; i<*n; i++){
	printf("b[%d] = ", i);
	scanf ("%lf", &Matrix[i][*n]);
    }
    return Matrix;
}




double * Insert_X(int n) {
    //data of matrix X are inserted by the user
    double *Matrix; 
    int i;
    Matrix = (double*)malloc(n*sizeof(double));
    printf("\nParakalw eisagete to dianysma X (p.x (1,1,1,1,1)^T gia enan 5x5 pentadiagwnio pinaka) : \n");
    for (i=0; i<n; i++) {
	printf("X[%d] = ", i);
	scanf("%lf", &Matrix[i]);
    }
    return Matrix;
}



int GivenMatrix() {
    //menu for the 2nd selection(specific matrices) of the main manu
    int number;
    do{
	printf("\n\nSpecific matrices LIST		\nEpilekse sygkekrimeno pinaka (apo ta dedomena ths ekfwnhshs):\n\t1.Pentadiagwnios pinakas me diastaseis 5X5 (adjustment 1.1)\n\t2. Pentadiagwnios pinakas me diastaseis 10X10 (adjustment 1.ii)\n\t3. Pentadiagwnios pinakas me distaseis NXN, opoy N, a, b, c kai d dinontai apo esas me times 100, 1000, 10000(for N) kai 0.1-1.9(for each of a,b,c,d)\n\n\t4. Epistrofh sto menu.\n");	
	
	printf("INPUT : ");
	scanf("%d", &number);
    }while((number<1)||(number>4));
    return number;
}



double **Matrix_5x5(int *n, double **X) {
    //returns a specific matrix (adjustment 1.i)
    double **A;
    *n=5;
    A = dispArray(*n, *n+1);
    //the one more column corresponds to the b vector
    A[0][0]=4; 	     A[0][1]=-0.3; 	A[0][2]=-0.4; 	 A[0][3]=0; 	   A[0][4]=0; 	     A[0][5]=3.3;
    A[1][0]=-0.2;    A[1][1]=4; 	A[1][2]=-0.3; 	 A[1][3]=-0.4; 	   A[1][4]=0; 	     A[1][5]=3.1;
    A[2][0]=-0.1;    A[2][1]=-0.2; 	A[2][2]=4; 	 A[2][3]=-0.3; 	   A[2][4]=-0.4;     A[2][5]=3;
    A[3][0]=0; 	     A[3][1]=-0.1; 	A[3][2]=-0.2; 	 A[3][3]=4; 	   A[3][4]=-0.3;     A[3][5]=3.4;
    A[4][0]=0;	     A[4][1]=0; 	A[4][2]=-0.1; 	 A[4][3]=-0.2; 	   A[4][4]=4; 	     A[4][5]=3.7;


    *X = (double*)malloc(*n*sizeof(double));
    (*X)[0]=1; 	   (*X)[1]=1; 	 (*X)[2]=1; 	(*X)[3]=1;     (*X)[4]=1;
    return A;
}

double **Matrix_10x10(int *n, double **X) {

    //returns a specific matrix (adjustment 1.ii)
    double **A;
    *n=10;

    A = dispArray(*n, *n+1);
    //the one more column corresponds to the b vector

    A[0][0]=4; 		A[0][1]=-0.2; 		A[0] [2]=-0.1; 	     A[0][3]=0; 	A[0][4]=0; 	     A[0] [5]=0; 	A[0][6]=0; 		A[0][7]=0; 	   A[0][8]=0; 		A[0][9]=0; 	  A[0][10]=3.7;
    A[1][0]=-0.3; 	A[1][1]=4; 		A[1][2]=-0.2; 	     A[1][3]=-0.1;	A[1][4]=0; 	     A[1][5]=0; 	A[1][6]=0; 		A[1][7]=0; 	   A[1][8]=0; 		A[1][9]=0;	  A[1][10]=3.4;
    A[2][0]=-0.4; 	A[2][1]=-0.3; 		A[2][2]=4; 	     A[2][3]=-0.2; 	A[2][4]=-0.1; 	     A[2][5]=0; 	A[2][6]=0; 		A[2][7]=0; 	   A[2][8]=0; 		A[2][9]=0;	  A[2][10]=3;
    A[3][0]=0; 		A[3][1]=-0.4; 		A[3][2]=-0.3; 	     A[3][3]=4; 	A[3][4]=-0.2; 	     A[3][5]=-0.1; 	A[3] [6]=0; 		A[3][7]=0;	   A[3][8]=0; 		A[3][9]=0; 	  A[3][10]=3;
    A[4][0]=0; 		A[4][1]=0; 		A[4][2]=-0.4; 	     A[4][3]=-0.3; 	A[4][4]=4; 	     A[4][5]=-0.2; 	A[4][6]=-0.1; 		A[4][7]=0; 	   A[4][8]=0;	    	A[4][9]=0; 	  A[4][10]=3;
    A[5][0]=0; 		A[5][1]=0; 		A[5][2]=0; 	     A[5][3]=-0.4; 	A[5] [4]=-0.3; 	     A[5][5]=4; 	A[5][6]=-0.2;		A[5][7]=0;	   A[5][8]=0; 		A[5][9]=0; 	  A[5][10]=3;	
    A[6][0]=0;		A[6][1]=0; 		A[6] [2]=0; 	     A[6][3]=0; 	A[6][4]=-0.4; 	     A[6] [5]=-0.3; 	A[6][6]=4; 		A[6][7]=-0.2;      A[6][8]=-0.1; 	A[6][9]=0;	  A[6][10]=3;					
    A[7][0]=0; 		A[7][1]=0;		A[7][2]=0; 	     A[7][3]=0;		A[7][4]=0; 	     A[7][5]=-0.4; 	A[7] [6]=-0.3; 	        A[7][7]=4; 	   A[7][8]=-0.2; 	A[7][9]=-0.1; 	  A[7][10]=3;		
    A[8][0]=0;		A[8][1]=0;		A[8] [2]=0;	     A[8][3]=0;		A[8][4]=0; 	     A[8][5]=0; 	A[8][6]=-0.4; 		A[8][7]=-0.3; 	   A[8][8]=4; 		A[8][9]=-0.2; 	  A[8][10]=3.1;
    A[9][0]=0; 		A[9][1]=0; 		A[9][2]=0; 	     A[9][3]=0; 	A[9][4]=0; 	     A[9][5]=0; 	A[9][6]=0;		A[9][7]=-0.4; 	   A[9][8]=-0.3; 	A[9][9]=4; 	  A[9][10]=3.3;


    *X = (double*)malloc(*n*sizeof(double));
    (*X)[0]=1; 		(*X)[1]=1; 	 (*X) [2]=1; 	(*X) [3]=1; 	(*X)[4]=1; 	(*X) [5]=1; 	(*X) [6]=1; 	(*X) [7]=1; 	(*X) [8]=1; 	(*X) [9]=1;
    return A;
}





double **Matrix_NxN (int *n, double **X) {
    //returns a specific matrix (adjustment 2)
    double **A, a, b, c, d;
    //where -a=k, -b=l, -c=r, -d=s

    int i, j;
    printf("Parakalw eisagete th diastash toy pentadiagwnioy pianaka[NXN]  (p.x. 100, 1000, 10000) : ");
    scanf("%d", n);
    printf("Parakalw eisagete tiw parametroys a, b, c, d antistoixa se aythn th seira me a, b, c, d sto diasthma 0.1(0.1)1.9 (p.x. a=1.2, b=0.9, c=0.6, d=0.3) :\n");
    printf("a = ");
    scanf("%lf", &a);
    printf(" b = ");
    scanf("%lf", &b);
    printf("c=	"); 
    scanf("%lf", &c);
    printf("d = ");
    scanf("%lf", &d);

    *X = (double*) malloc(*n*sizeof(double));
    A = dispArray(*n, *n+1);
    //the one more column corresponds to the b vector
    //values for matrix A
    for(i=0; i<*n; i++)
	for(j=0; j<*n; j++)
	    switch(abs(i-j)){
		case 0: A[i][j]=4;
			break;
		case 1: if(i>j) A[i][j]=-b; 
			else A[i][j]=-c;
			break;
		case 2: if (i>j) A[i][j]=-a;
			else A[i][j]=-d;
			break;
		default:A[i][j]=0.0;
			break;
	    }

    //values for x

    for(i=0; i<*n; i++)
	(*X)[i]=1;

    //b computation
    for(i=0; i<*n; i++) {
	A[i][*n] = 0;  //initializes into zero all the values of b vector
	for(j=0; j<*n; j++)
	    A[i][*n] += A[i][j]*(*X) [j];
	    //the above vector b cause of x=(1,1, ...,1)^T means b[i] = a[i]+b[i]+diag[i]+c[i]+d[i] = -(k[i]+l[i]+d[i]+r[i]+s[i])
    }
    return A;
}


double ** Insert_rand_A(int *n) {
    //data of matrix A are randomly selected and inserted (adjustment 3)
    srand(time(NULL));	//where	k=-a, l=-b, r=-c, S=-d
    double **A, k, l, r, s, values, diag_value, mult;
    int i, j, case_values;

    printf("Give the size N of the pentadiagonial matrix[NXN] (e.g. 50, 100, 1000) :  ");
    scanf("%d",n);
    //values for parameters of matrix A
    values = 0.1+(0.1* (double)(rand()%19));		//desirable values of parameters a,b,c,d between [0.1,1.9]
    case_values = (int) (values* 10);
    diag_value = (sign()*(4.0+(double)((rand()%385)/4)))-1;		//desirable values of diagonial between [-100,-4] or [4,100] all multiplies of 4(-100,...,-4,4,8,12, ...,


    switch(case_values) {
        case 1:
        case 2:
	case 3:
	case 4:
	case 5:
		mult = 0.1*(1.0+(double)(rand()%4));
		k = sign() *values;
		l = sign()* (k+mult);
		r = sign()*(k+2.0*mult);
		s = sign()*(k+3.0*mult);
		break;
	case 6:
	case 7:
	case 8:
	case 9:
	case 10:
		mult = 0.1*(1.0+(double)(rand()%3));
		k = sign()*values;
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
    printf("\n\n\nPentadiagonial matrix parameters selected as follows:\nk=%2.1lf, l=%2.1lf, d=%2.1lf, r=%2.11f, s=%2.1lf\n\n\n", k, l, diag_value,r,s);
    A = dispArray(*n, *n+1);
    //the one more column corresponds to the b vector
    //values for matrix A

    for(i=0; i<*n; i++)
	for(j=0; j<*n; j++)
	    switch(abs(i-j)) {
		case 0: A[i][j]=diag_value;
			break;
		case 1: if (i>j) A[i][j]=l ;
			else A[i][j]=r;
			break;
		case 2: if (i>j) A[i][j]=k;
			else A[i][j]=s;
			break;
		default:A[i][j]=0.0;
			break;
		}
    //randomly selected data for b vector below so that X is (1,1, ...,1,1)^T
    for(i=0; i<*n; i++) {
	A[i][*n]=0.0;	//initializes into all the values of b vector
        for(j=0; j<*n; j++)
	    A[i][*n] += A[i][j];
        //the above vector b cause of x=(1,1, ...,1)^T means b[i] = a[i]+b[i]+diag[i]+[i]+d[i] = -(k[i]+[i]+d[i]+r[i]+s[i])
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

    if((file=fopen("ask2_PSD.txt", "r"))==NULL) {
        perror("FILE ERROR : Either the file with name 'ask2_PSD.txt' does not exist or it is placed in a wrong directory");
        exit(0);
        //not return 0..
    }

    fscanf(file, "%d", n);	//the first line in file includes only the number of A dimensions
    printf("Dimensions are %dx%d\n", *n, *n);
    A = dispArray(*n, *n+1);
    //the one more column corresponds to the b vector
    (*X) = (double*)malloc(*n*sizeof(double));

    //reading data from file and import them in matrix A below
    for(i=0; i<*n; i++)
	for(j=0; j<*n; j++)
	    fscanf(file, "%lf",&A[i][j]);

    //reading data from file and import them in vector b below
    for (i=0; i<*n; i++)
	fscanf(file, "%lf", &A[i][*n]);

    //reading data from file and import them in vector x below
    for(i=0; i<*n; i++)
	fscanf(file, "%lf",&(*X)[i]);

    fclose(file);
    return A;
}




double ** dispArray(int n, int m) {
    double **Matrix;
    int i, j;

    Matrix= (double**)malloc(n*sizeof(double*));

    for(i=0;i<n;i++)
        *(Matrix+i)=(double*) malloc(m*sizeof(double));
    return Matrix;
}


double sign() {
    if(rand()%2==0)
        return 1;
    else
    	return -1.;
}

