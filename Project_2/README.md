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




<br> <br> <br>



<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>
<h1>Part 2</h1>

<p>This program is designed to solve a linear system of equations ( Ax = b ) using the Partially Separated Decomposition (PSD) method. The program provides various options for inputting the matrix ( A ) and vector ( b ), including manual entry, selection from predefined matrices, random generation within specified parameters, or reading from a file. After inputting the data, the program allows for the execution of the PSD method with options to experimentally adjust the accuracy or study the convergence of the iterative method with different parameters.</p>


<h2>Compilation and Execution:</h2>

<p>To compile the code, use a C compiler such as gcc:</p>

<pre><code>gcc -o particle particle.c -lm</code></pre>

<p>To execute the compiled program:</p>

<pre><code>./particle</code></pre>



<h2>Usage</h2>

<h3>Compilation</h3>

<p>To compile the code, use any standard C compiler. For example, you can use <code>gcc</code>:</p>

<h3>Execution</h3>

<p>After compilation, run the executable <code>main</code>. The program will display a menu with options for inputting the matrix ( A ) and vector ( b ) and choosing the method of execution.</p>

<h2>Functionality</h2>

<h3>Input Methods</h3>

<ol>
<li><strong>Manual Input</strong>: Allows the user to manually input the dimensions and elements of matrix \( A \) and vector \( b \).</li>
<li><strong>Predefined Matrices</strong>: Provides a list of predefined matrices, including a 5x5, 10x10, or user-defined \( N \times N \) matrix with specified parameters.</li>
<li><strong>Random Matrix</strong>: Generates a random \( N \times N \) matrix within specified parameter ranges.</li>
<li><strong>File Input</strong>: Reads matrix \( A \) and vector \( b \) from a file named 'ask2_ESOR.txt'.</li>
</ol>

<h3>PSD Method Execution</h3>

<ol>
<li><strong>Experimental Adjustment</strong>: Executes the PSD method with specified parameters \( t \) and \( w \), allowing the user to adjust the accuracy.</li>
<li><strong>Convergence Study</strong>: Conducts an experimental study of the convergence of the PSD iterative method with parameters \( t \) and \( w \) ranging from 0.1 to 1.9 with a step of 0.1.</li>
<li><strong>Exit</strong>: Terminates the program.</li>
</ol>

<h2>File Structure</h2>

<ul>
<li><code>main.c</code>: Contains the main code implementing the PSD method and menu functionalities.</li>
<li><code>ask2_ESOR.txt</code>: Sample file for reading matrix \( A \) and vector \( b \).</li>
<li><code>README.md</code>: This file, providing instructions and information about the code.</li>
</ul>

<h2>Functions:</h2>

<ul>
    <li><code>menu()</code>:
        <ul>
            <li>This function displays the main menu for the program and prompts the user to select an option.</li>
        </ul>
    </li>
    <li><code>Insert_A(int*)</code>:
        <ul>
            <li>Prompts the user to input the dimension N of the matrix A and then inserts the elements of the matrix along with the elements of the vector b.</li>
        </ul>
    </li>
    <li><code>Insert_X(int)</code>:
        <ul>
            <li>Prompts the user to input the vector X (initial guess for the solution).</li>
        </ul>
    </li>
    <li><code>GivenMatrix()</code>:
        <ul>
            <li>Displays a menu for selecting specific pre-defined matrices or returning to the main menu.</li>
        </ul>
    </li>
    <li><code>Matrix_5x5(int*, double**)</code>, <code>Matrix_10x10(int*, double**)</code>, <code>Matrix_NxN(int*, double**)</code>:
        <ul>
            <li>Generate specific matrices (5x5, 10x10, or NxN) along with the corresponding vector X.</li>
        </ul>
    </li>
    <li><code>Read_file_Axb(int*, double**)</code>:
        <ul>
            <li>Reads a matrix A and the corresponding vector b from a file named 'ask2_ESOR.txt'.</li>
        </ul>
    </li>
    <li><code>PSD_menu()</code>:
        <ul>
            <li>Displays the menu for selecting the method of implementing the PSD algorithm.</li>
        </ul>
    </li>
    <li><code>PSD_mainmethod(double**, double*, int)</code>:
        <ul>
            <li>Implements the PSD method to solve the linear system Ax = b. It prompts the user to input parameters for the PSD method and executes the method accordingly.</li>
        </ul>
    </li>
    <li><code>PSD_iterativeMethod(double**, double**, int, double, double, int, int*)</code>:
        <ul>
            <li>Implements the iterative process of the PSD method for factorizing matrix A into L and U using partial pivoting.</li>
        </ul>
    </li>
    <li><code>dispArray(int, int)</code>:
        <ul>
            <li>Dynamically allocates memory for a 2D array of given dimensions.</li>
        </ul>
    </li>
</ul>



<h2>Notes</h2>
<ul>
<li>The program relies on standard C libraries and is compatible with most C compilers.</li>
<li>Ensure proper input format and follow the instructions provided by the program for accurate results.</li>
<li>Experiment with different input methods and parameters to understand the behavior of the PSD method in solving linear systems.</li>
</ul>

</body>
</html>
