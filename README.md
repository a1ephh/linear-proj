This is a bonus component implemented in c++, our code provides a diagnostic tool to verify if a given mathematical function satisfies the axioms of Inner Product Space. We can use it to check for positivity, symmetry and linearity across different vector spaces including vectors and polynomials 

Main features: 
Axiom validation: 
Positivity: <v,v>  >=0 
Symmetry: <u,v> == <v,u> 
Linearity: <au + bv,w> = a <u,w> + b <v,w>

Implementations of the following:
Standard Dot Product 
Weighted Inner Product
Polynomial Integrals 
Code structure: 
Main.cpp: Contains the logic to run the test cass
Implementation.h: Defines data structures for vectors, R^3 and polynomials
Utils.h: Contains the operators for vector and polynomial addition, subtraction and scalar multiplication

How to run:

This project requires a C++ compiler supporting the C++17 standard (necessary for structured bindings used in the linearity tests).

Run this in the terminal:
g++ -std=c++17 -o ipchecker main.cpp


The program accepts several terminal arguments to customize the test suite:

1. Standard Dot Product (Default n=3):
   Run:
   ./ipchecker --function dot

2. Weighted Inner Product (n=2 with custom weights):
   Run:
   ./ipchecker -d 2 -f weighted -w 2.0 5.0

3. Invalid "Fail" Cases (Testing non-examples):
   Run:
   ./ipchecker --function custom

4. Run All Tests:
   Run:
   ./ipchecker --function all

Understanding the Output
For every function tested, the program outputs:
Header: Displays the function type, dimension, and tolerance used.
Step-by-Step Results: Shows the calculated values for both sides of the axiom equations (e.g., <u,v> vs <v,u>) and the calculated difference.
Status: A clear PASS/FAIL label for each individual test case.  
Summary Report: An overall conclusion on whether the function is a valid inner product based on the collective axiom results. 
