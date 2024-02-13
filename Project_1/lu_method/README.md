<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>

<h1>Linear Algebra Solver with LU Decomposition in C</h1>

<p>This project consists of a collection of C functions for performing various linear algebra operations, including LU decomposition, matrix inversion, matrix multiplication, and more. The project also includes a main program demonstrating the usage of these functions to solve linear systems and compute errors in solutions.</p>

<h2>Contents:</h2>

<ul>
  <li><code>LU_substitution</code>: Performs forward and backward substitutions using LU-decomposed matrices to solve linear systems.</li>
  <li><code>LU_decompose</code>: Performs LU decomposition of a matrix using Crout's method with partial implicit pivoting.</li>
  <li><code>invert_matrix</code>: Uses LU decomposition and back substitution to invert a matrix.</li>
  <li><code>MatMult</code>: Multiplies two square real matrices.</li>
  <li><code>MatPrint</code>: Prints a square real matrix.</li>
  <li><code>MatMultOneDim</code>: Multiplies a square real matrix by a one-dimensional matrix.</li>
  <li><code>MatSubtraction</code>: Subtracts two real matrices.</li>
  <li><code>MatSubtraction2</code>: Subtracts two square real matrices.</li>
  <li><code>MatRowsSum</code>: Sums the rows of a square real matrix.</li>
  <li><code>readmatrix</code>: Reads matrix data from a file.</li>
  <li><code>areDoubleEqual</code>: Checks if two double values are equal.</li>
  <li><code>compareDoubles</code>: Compares two double values for sorting purposes.</li>
</ul>

  <h2>Prerequisites</h2>
  <ul>
    <li>C compiler (e.g., GCC)</li>
  </ul>

<h2>Usage:</h2>

<ol>
  <li>Clone the repository.</li>
  <li>Include the necessary header files in your C code.</li>
  <li>Utilize the provided functions in your program for various linear algebra operations.</li>
</ol>

<h2>Input:</h2>

<ul>
  <li>The input matrices are typically read from files using the <code>readmatrix</code> function.</li>
  <li>Alternatively, matrices can be initialized manually within the code.</li>
</ul>

<h2>Output:</h2>

<ul>
  <li>Solutions to linear systems are printed using <code>MatPrint</code>.</li>
  <li>Absolute errors in solutions are calculated and printed within the main program.</li>
</ul>



</body>
</html>






<br> <br> <br>
##**<p align=center>Advanced Matrix Operations with LU Decomposition in C**


<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>
  <h1>Efficient Matrix Operations with LU Decomposition for Inversion/h1>

  <h2>Overview</h2>
  <p>This project provides a collection of functions for performing various matrix operations, including LU decomposition, matrix inversion, matrix multiplication, and more. The functions are written in C and can be utilized in any C program that requires matrix manipulation.</p>

  <h2>Features</h2>
  <ul>
    <li>LU decomposition of a square matrix</li>
    <li>Matrix inversion using LU decomposition</li>
    <li>Matrix multiplication</li>
    <li>Subtraction of matrices</li>
    <li>Summation of rows in a matrix</li>
    <li>Reading matrices from files</li>
    <li>Error calculation for matrix inversion accuracy</li>
    <li>Utilizes Cholesky method for matrix inversion</li>
  </ul>

  <h2>Prerequisites</h2>
  <ul>
    <li>C compiler (e.g., GCC)</li>
  </ul>

  <h2>Usage</h2>
  <ol>
    <li>Clone the repository or download the source code files.</li>
    <li>Include the necessary header files (<code>stdio.h</code>, <code>stdlib.h</code>, <code>math.h</code>, <code>stddef.h</code>, <code>float.h</code>) in your C program.</li>
    <li>Include the provided header file <code>matrix_operations.h</code> in your C program to access the functions.</li>
    <li>Compile your program along with the source files <code>matrix_operations.c</code>.</li>
    <li>Utilize the provided functions in your C program for matrix manipulation as needed.</li>
  </ol>

  <h2>Functionality</h2>
  <ul>
    <li><strong>LU Decomposition:</strong> Decompose a square matrix into lower and upper triangular matrices using Crout's method with partial implicit pivoting.</li>
    <li><strong>Matrix Inversion:</strong> Invert a square matrix using LU decomposition and back substitution. Supports error checking for singular matrices. Calculates the absolute relative error for inversion accuracy.</li>
    <li><strong>Matrix Multiplication:</strong> Multiply two square matrices to produce a resulting matrix.</li>
    <li><strong>Matrix Subtraction:</strong> Subtract one matrix from another to produce a resulting matrix.</li>
    <li><strong>Summation of Rows:</strong> Calculate the sum of each row in a square matrix.</li>
    <li><strong>File Input:</strong> Read matrices from files for processing.</li>
  </ul>

  <h2>Input</h2>
  <ul>
    <li>Matrices can be read from input files in a specific format (e.g., space-separated values).</li>
    <li>The size and content of matrices should adhere to the defined constraints.</li>
  </ul>

  <h2>Output</h2>
  <ul>
    <li>Functions produce output matrices or numerical results based on the operation performed.</li>
    <li>Output may include inverse matrices, multiplied matrices, error calculations, or other relevant results.</li>
  </ul>

  <h2>Error Calculation</h2>
  <ul>
    <li>The library includes functionality to calculate the absolute relative error for matrix inversion accuracy.</li>
    <li>Error calculation is performed by comparing the original matrix with the inverted matrix.</li>
  </ul>

</body>
</html>

