<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>

<p align=center> <h1>1.Linear System Solver with Cholesky Method</h1>

<p>This part implements a linear system solver using the Cholesky decomposition method. The Cholesky decomposition is particularly useful for solving linear systems with symmetric positive definite matrices. The program takes a matrix A and a vector B as input and solves the equation Ax = B, where A is symmetric positive definite.</p>

<h2>Overview</h2>

<p>The project consists of several functions:</p>

<ul>
  <li><code>choldc</code>: Performs Cholesky decomposition on a given matrix A.</li>
  <li><code>choldcsl</code>: Computes the inverse of a lower triangular Cholesky decomposed matrix.</li>
  <li><code>choldet</code>: Computes the determinant of a symmetric positive definite matrix using Cholesky decomposition.</li>
  <li><code>cholsl</code>: Computes the inverse of a symmetric positive definite matrix using Cholesky decomposition.</li>
  <li><code>func</code>: Solves the linear system Ax = B using Gaussian elimination method.</li>
  <li><code>LU_substitution</code>: Performs forward and backward substitution using LU decomposition.</li>
  <li>Other utility functions for matrix operations such as multiplication, subtraction, printing, and checking positive definiteness.</li>
</ul>

<h2>Usage</h2>

<ol>
  <li><strong>Compilation:</strong> Compile the program using a C compiler such as <code>gcc</code>.</li>
  
  <pre><code>gcc -o cholesky_linear_solving cholesky_linear_solving.c -lm</code></pre>

  <li><strong>Execution:</strong> Run the compiled program with the required input.</li>
  
  <pre><code>./cholesky_linear_solving</code></pre>

  <li><strong>Input:</strong> The program reads a matrix A and a vector B from a file named <code>2a1.txt</code>. Make sure this file exists in the same directory as the executable.</li>

  <li><strong>Output:</strong> The program prints the original matrix A, solves the linear system using the Cholesky method, and displays the results including the L and U matrices, the solution vector x, and the absolute error and remainder.</li>
</ol>

<h2>File Structure</h2>

<ul>
  <li><code>cholesky_linear_solving.c</code>: Main source code file containing the implementation of the linear system solver.</li>
  <li><code>2a1.txt</code>: Input file containing the matrix A for the linear system.</li>
</ul>

<h2>Additional Notes</h2>

<p>The program assumes that the input matrix A is symmetric positive definite. It may not produce correct results if the input matrix does not meet this condition.</p>






<br> <br> <br>

<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>

<p align=center> <h1>2.Matrix Inversion using Cholesky Decomposition</h1>

<p>This part implements the inversion of a square real symmetric matrix using the Cholesky decomposition method. The algorithm is designed for matrices that are positive definite.</p>

<h2>Description</h2>

<p>The Cholesky decomposition method is employed to find the inverse of the input matrix. The process involves decomposing the original matrix into a lower triangular matrix and its transpose, and then using these decomposed matrices to find the inverse. The resulting inverse matrix is used to verify the correctness of the inversion by multiplying it with the original matrix to obtain the identity matrix.</p>

<h2>Functions</h2>

<h3><code>cholesky_decomposition()</code></h3>

<p>This function performs the Cholesky decomposition of a symmetric positive definite matrix. It takes the following parameters:</p>

<ul>
  <li><code>double A[][]</code>: The input matrix to be decomposed.</li>
  <li><code>int n</code>: The size of the input matrix.</li>
  <li><code>double L[][]</code>: The resulting lower triangular matrix L.</li>
</ul>

<p>The function modifies the lower triangular matrix <code>L</code> in place.</p>

<!-- Other functions omitted for brevity -->
<h3><code>cholesky_decomposition()</code></h3>

<p>This function performs the Cholesky decomposition of a symmetric positive definite matrix. It takes the following parameters:</p>

<ul>
  <li><code>double A[][]</code>: The input matrix to be decomposed.</li>
  <li><code>int n</code>: The size of the input matrix.</li>
  <li><code>double L[][]</code>: The resulting lower triangular matrix L.</li>
</ul>

<p>The function modifies the lower triangular matrix <code>L</code> in place.</p>

<h3><code>forward_substitution()</code></h3>

<p>This function performs forward substitution to solve a lower triangular system of equations. It takes the following parameters:</p>

<ul>
  <li><code>double L[][]</code>: The lower triangular matrix from the Cholesky decomposition.</li>
  <li><code>double b[]</code>: The right-hand side vector of the system of equations.</li>
  <li><code>double x[]</code>: The solution vector to be computed.</li>
</ul>

<p>The function modifies the solution vector <code>x</code> in place.</p>

<h3><code>backward_substitution()</code></h3>

<p>This function performs backward substitution to solve an upper triangular system of equations. It takes the following parameters:</p>

<ul>
  <li><code>double U[][]</code>: The upper triangular matrix obtained from the Cholesky decomposition.</li>
  <li><code>double y[]</code>: The right-hand side vector of the system of equations (obtained from forward substitution).</li>
  <li><code>double x[]</code>: The solution vector to be computed.</li>
</ul>

<p>The function modifies the solution vector <code>x</code> in place.</p>


<h2>Usage</h2>

<p>To use this program, follow these steps:</p>

<ol>
  <li>Compile the program using a C compiler (e.g., GCC).</li>
  <pre><code>gcc cholesky_reversed.c -o cholesky_reversed -lm</code></pre>
  <li>Run the executable with the appropriate inputs.</li>
  <li>Provide the size of the matrix and the matrix itself.</li>
  <li>The program will output the original matrix, its inverse, and the result of multiplying the original matrix with its inverse.</li>
</ol>

<h2>Sample Output</h2>

<pre><code>
Inversion of a square real symmetric matrix by Cholesky method
(The matrix must be positive definite).

Size of the matrix is 4:   
The MATRIX is:

        |   1   1   1   1   |
        |   1   2   3   4   |
        |   1   3   6   10  |
        |   1   4   10  20  |

Find A^(-1) with Cholesky method

Matrix A:
1.000000   1.000000   1.000000   1.000000   
1.000000   2.000000   3.000000   4.000000   
1.000000   3.000000   6.000000   10.000000   
1.000000   4.000000   10.000000  20.000000   

Matrix Inv(A):
8.000000   -28.000000  56.000000  -70.000000  
-28.000000  140.000000 -322.000000 434.000000  
56.000000  -322.000000 812.000000 -1162.000000 
-70.000000 434.000000  -1162.000000 1742.000000 

Here is a verification:
Verification A * Inv(A) = I:
1.000000   0.000000   0.000000   0.000000   
0.000000   1.000000   0.000000   0.000000   
0.000000   0.000000   1.000000   0.000000   
0.000000   0.000000   0.000000   1.000000 
</code></pre>

<h2>Error Calculation</h2>

<p>The program also calculates the absolute relative error between the original matrix and its inverse, as well as the error resulting from multiplying the original matrix with its inverse and comparing the result with the identity matrix.</p>

<h2>File Input</h2>

<p>The program supports reading the input matrix from a file. Provide the filename containing the matrix values in the format specified by the program.</p>

</body>
</html>

</body>
</html>
