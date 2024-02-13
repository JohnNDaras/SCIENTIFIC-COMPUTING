<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Linear System Solver with Cholesky Method</title>
</head>
<body>

<h1>Linear System Solver with Cholesky Method</h1>

<p>This project implements a linear system solver using the Cholesky decomposition method. The Cholesky decomposition is particularly useful for solving linear systems with symmetric positive definite matrices. The program takes a matrix A and a vector B as input and solves the equation Ax = B, where A is symmetric positive definite.</p>

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
  
  <pre><code>gcc -o linear_solver linear_solver.c -lm</code></pre>

  <li><strong>Execution:</strong> Run the compiled program with the required input.</li>
  
  <pre><code>./linear_solver</code></pre>

  <li><strong>Input:</strong> The program reads a matrix A and a vector B from a file named <code>2a1.txt</code>. Make sure this file exists in the same directory as the executable.</li>

  <li><strong>Output:</strong> The program prints the original matrix A, solves the linear system using the Cholesky method, and displays the results including the L and U matrices, the solution vector x, and the absolute error and remainder.</li>
</ol>

<h2>File Structure</h2>

<ul>
  <li><code>linear_solver.c</code>: Main source code file containing the implementation of the linear system solver.</li>
  <li><code>2a1.txt</code>: Input file containing the matrix A for the linear system.</li>
</ul>

<h2>Additional Notes</h2>

<p>The program assumes that the input matrix A is symmetric positive definite. It may not produce correct results if the input matrix does not meet this condition.</p>

</body>
</html>
