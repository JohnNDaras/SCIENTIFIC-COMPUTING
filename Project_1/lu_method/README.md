<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Linear Algebra Operations in C</title>
</head>
<body>

<h1>Linear Algebra Operations in C</h1>

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
