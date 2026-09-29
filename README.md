# Matrix Calculator with Advanced Operations

A simple C++ console application for performing common matrix operations, including basic arithmetic and linear algebra utilities.

## Features

This project supports the following operations:

- Matrix addition
- Matrix subtraction
- Matrix multiplication
- Transpose
- Determinant
- Inverse
- Row-reduced echelon form (RREF)
- Rank

## Built With

- C++17
- Standard Template Library (STL)
- Console-based input/output

## How It Works

The program runs in a menu-driven format.  
You enter the dimensions and values of matrices, then choose an operation from the menu.

It handles invalid inputs and checks for cases such as:
- dimension mismatch
- non-square matrices for determinant/inverse
- singular matrices with no inverse

## Operations

### 1. Addition
Adds two matrices of the same size.

### 2. Subtraction
Subtracts one matrix from another of the same size.

### 3. Multiplication
Multiplies matrix `A` by matrix `B` when the number of columns in `A` matches the number of rows in `B`.

### 4. Transpose
Reverses rows and columns of a matrix.

### 5. Determinant
Computes the determinant of a square matrix using Gaussian elimination with partial pivoting.

### 6. Inverse
Computes the inverse of a square matrix using Gauss-Jordan elimination.  
If the matrix is singular, the program reports that no inverse exists.

### 7. RREF
Converts the matrix into reduced row echelon form.

### 8. Rank
Finds the rank of the matrix using RREF.

## How to Compile

Use the following command:

```bash
g++ -std=c++17 -O2 -o matrix_calculator matrix_calculator.cpp
```

## How to Run

After compiling, run:

```bash
./matrix_calculator
```

## Example Use Case

This project is useful for:
- practicing C++
- learning matrix operations
- exploring linear algebra concepts
- building a beginner-friendly mathematics tool

## Project Structure

```text
matrix_calculator.cpp
README.md
```

## Notes

- This is a console-based project.
- All calculations are done using `double`.
- A small epsilon value is used to handle floating-point precision issues.
- The program prints matrices in a formatted way for better readability.

## Future Improvements

Possible extensions for this project:
- support for matrix input from files
- determinant by cofactor expansion for small matrices
- eigenvalue/eigenvector tools
- solving linear systems
- exporting results

## Author

Made by Arpit
