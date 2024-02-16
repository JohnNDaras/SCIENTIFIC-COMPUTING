<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>

<h1>Part 1</h1>

<p>This program implements the ESOR (Elementary Symmetric Overrelaxation) iterative method for solving linear systems of equations of the form Ax = b. It provides various options for inputting the matrix A and vector b, including manual input, predefined matrices, random generation, and reading from a file. Additionally, it offers an experimental validation of the method's accuracy and a study of the convergence behavior by varying the parameters within a specified range.</p>

<h2>Compilation and Execution:</h2>

<p>To compile the code, use a C compiler such as gcc:</p>

<pre><code>gcc -o extrapolated extrapolated.c -lm</code></pre>

<p>To execute the compiled program:</p>

<pre><code>./extrapolated</code></pre>

<h2>Usage:</h2>

<p>Upon execution, the program will display a menu allowing you to choose among several options:</p>

<ul>
    <li>Manual Input: Enter the matrix A and vector b manually.</li>
    <li>Predefined Matrices: Choose from a list of predefined matrices.</li>
    <li>Random Matrix: Generate a random matrix A and vector b.</li>
    <li>Read from File: Read matrix A and vector b from a file named 'ask2_ESOR.txt'.</li>
    <li>Exit: Terminate the program.</li>
</ul>

<h3>File Format for Input:</h3>

<p>If you choose to read from a file, ensure the file 'ask2_ESOR.txt' follows this format:</p>

<ul>
    <li>First line: Integer representing the dimension N of the matrix A (N x N).</li>
    <li>Next N lines: Matrix A values separated by spaces.</li>
    <li>Next line: Vector b values separated by spaces.</li>
    <li>Final line: Vector X values separated by spaces.</li>
</ul>

<p>Note:</p>

<ul>
    <li>The ESOR algorithm requires specifying parameters t and w, which can be entered manually or determined experimentally for optimal convergence.</li>
    <li>The program provides feedback on the execution time and the solution of the linear system.</li>
</ul>

<h2>Main Functions:</h2>
<ul>
        <li><code>main()</code></li>
        <p>The main function presents a menu-driven interface for selecting the method of input for the coefficient matrix A and vector b. It then calls the <code>ESOR_mainmethod()</code> function to perform the ESOR iterative method for solving the linear system.</p>
        <li><code>ESOR_mainmethod()</code></li>
        <p>This function is the main driver for performing the ESOR iterative method. It takes the coefficient matrix A, the vector of unknowns x, and the size of the matrix as input. It also prompts the user for the maximum number of iterations and provides options for experimental validation and convergence study of the ESOR method.</p>
        <li><code>ESOR_iterativeMethod()</code></li>
        <p>This function implements the ESOR iterative method to solve the linear system Ax = b. It takes the coefficient matrix A, an initial approximation of the solution x1, the size of the matrix, relaxation parameters (t and w), maximum number of iterations, and a pointer to store the iteration count as input.</p>
    </ul>
    <h2>Matrix Generation Functions:</h2>
    <p>These functions (<code>Insert_A()</code>, <code>Insert_X()</code>, <code>GivenMatrix()</code>, <code>Matrix_5x5()</code>, <code>Matrix_10x10()</code>, <code>Matrix_NxN()</code>, <code>Insert_rand_A()</code>, <code>Insert_rand_X()</code>, <code>Read_file_Axb()</code>) are responsible for generating various types of coefficient matrices A and vectors b. They provide options for manual input, selection from predefined matrices, random generation, and reading from a file.</p>
    <h2>Utility Functions:</h2>
    <ul>
        <li><code>dispArray()</code>: Allocates memory for a 2D array (matrix).</li>
        <li><code>sign()</code>: Generates a random sign (+1 or -1).</li>
    </ul>


</body>
</html>
