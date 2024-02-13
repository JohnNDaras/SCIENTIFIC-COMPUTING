<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Linear Algebra Operations in C</title>
</head>
<body>

<h1>Linear Algebra Operations in C</h1>

<p>This C project provides functionalities for various linear algebra operations including LU decomposition, matrix inversion, matrix multiplication, and more. The project is structured into several functions each serving a specific purpose.</p>

<h2>Contents</h2>

<ol>
<li><a href="#introduction">Introduction</a></li>
<li><a href="#functions">Functions</a></li>
<li><a href="#usage">Usage</a></li>
<li><a href="#compilation">Compilation</a></li>
<li><a href="#examples">Examples</a></li>
</ol>

<h2 id="introduction">Introduction</h2>

<p>The project includes functions for solving linear systems of equations, matrix inversion, matrix multiplication, and other related operations. It's particularly useful for applications in numerical analysis, computational mathematics, and scientific computing.</p>

<h2 id="functions">Functions</h2>

<ul>
<li><code>LU_decompose</code>: Performs LU decomposition of a matrix using Crout's method with partial implicit pivoting.</li>
<li><code>LU_substitution</code>: Performs forward and backward substitutions using LU-decomposed matrices to solve linear systems.</li>
<li><code>invert_matrix</code>: Utilizes LU decomposition and back substitution to invert a square matrix.</li>
<li><code>MatMult</code>: Multiplies two square matrices.</li>
<li><code>MatMultOneDim</code>: Multiplies a square matrix by a one-dimensional matrix.</li>
<li><code>MatSubtraction</code>: Subtracts two matrices.</li>
<li><code>MatSubtraction2</code>: Subtracts two square matrices.</li>
<li><code>MatRowsSum</code>: Calculates the sum of rows of a square matrix.</li>
<li><code>MatPrint</code>: Prints the contents of a matrix.</li>
<li><code>readmatrix</code>: Reads a matrix from a file.</li>
<li><code>areDoubleEqual</code>: Checks if two double precision floating-point numbers are equal.</li>
<li><code>compareDoubles</code>: Compares two double precision floating-point numbers for sorting purposes.</li>
</ul>

<h2 id="usage">Usage</h2>

<p>To use these functions, include the necessary headers and call the respective functions with appropriate arguments. Ensure that the matrix dimensions are compatible and handle errors gracefully.</p>

<h2 id="compilation">Compilation</h2>

<p>Compile the project using a C compiler such as GCC:</p>

<pre><code>gcc main.c -o linear_algebra_operations -lm
</code></pre>

<p>This command compiles the main file along with necessary libraries and generates an executable named <code>linear_algebra_operations</code>.</p>

<h2 id="examples">Examples</h2>

<p>Below is an example demonstrating how to use the provided functions:</p>

<pre><code>&lt;?php
#include &lt;stdio.h>
#include &lt;stdlib.h>
#include "linear_algebra_operations.h"

int main() {
    // Example usage of functions
    double A[NDIM][NDIM], B[NDIM] = {25.0, 75.0, 37.0, 46.0, 5.0, 93.0, -16.0, 41.0};
    // Initialize matrix A
    // Perform LU decomposition
    // Solve linear system Ax = B
    // Print solution

    return 0;
}
</code></pre>

<p>Note: Ensure that the NDIM macro is defined appropriately to match the dimensions of your matrices.</p>

</body>
</html>
