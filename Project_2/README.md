<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Linear System Resolution using ESOR Iterative Method</title>
</head>
<body>

<h1>Linear System Resolution using ESOR Iterative Method</h1>

<p>This program implements the ESOR (Elementary Symmetric Overrelaxation) iterative method for solving linear systems of equations of the form Ax = b. It provides various options for inputting the matrix A and vector b, including manual input, predefined matrices, random generation, and reading from a file. Additionally, it offers an experimental validation of the method's accuracy and a study of the convergence behavior by varying the parameters within a specified range.</p>

<h2>Compilation and Execution:</h2>

<p>To compile the code, use a C compiler such as gcc:</p>

<pre><code>gcc -o esor_linear_system_resolution esor_linear_system_resolution.c -lm</code></pre>

<p>To execute the compiled program:</p>

<pre><code>./esor_linear_system_resolution</code></pre>

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

<h2>License:</h2>

<p>This program is provided under the MIT License. Feel free to modify and distribute it according to your needs.</p>

</body>
</html>
