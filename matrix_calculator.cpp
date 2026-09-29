// Matrix Calculator with Advanced Operations
// Features: add, subtract, multiply, transpose, determinant, inverse,
//           row-reduction (RREF), rank
// Compile: g++ -std=c++17 -O2 -o matrix_calculator matrix_calculator.cpp

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>
using namespace std;

const double EPS = 1e-9;

class Matrix {
public:
    int rows = 0, cols = 0;
    vector<vector<double>> a;

    Matrix() = default;
    Matrix(int r, int c, double v = 0.0) : rows(r), cols(c), a(r, vector<double>(c, v)) {}

    static Matrix identity(int n) {
        Matrix I(n, n);
        for (int i = 0; i < n; i++) I.a[i][i] = 1.0;
        return I;
    }

    void read(const string& name) {
        cout << "Enter dimensions of " << name << " (rows cols): ";
        while (!(cin >> rows >> cols) || rows <= 0 || cols <= 0) {
            cin.clear();
            cin.ignore(10000, '\n');
            cout << "Invalid. Enter positive rows and cols: ";
        }
        a.assign(rows, vector<double>(cols));
        cout << "Enter " << rows << " x " << cols << " entries, row by row:\n";
        for (auto& row : a)
            for (double& x : row) cin >> x;
    }

    void print(const string& title = "") const {
        if (!title.empty()) cout << title << "\n";
        for (const auto& row : a) {
            cout << "  [ ";
            for (double x : row) {
                if (fabs(x) < EPS) x = 0.0;  // avoid printing -0
                cout << setw(10) << fixed << setprecision(4) << x << " ";
            }
            cout << "]\n";
        }
    }

    Matrix operator+(const Matrix& o) const {
        if (rows != o.rows || cols != o.cols) throw runtime_error("Dimension mismatch for addition.");
        Matrix r(rows, cols);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) r.a[i][j] = a[i][j] + o.a[i][j];
        return r;
    }

    Matrix operator-(const Matrix& o) const {
        if (rows != o.rows || cols != o.cols) throw runtime_error("Dimension mismatch for subtraction.");
        Matrix r(rows, cols);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) r.a[i][j] = a[i][j] - o.a[i][j];
        return r;
    }

    Matrix operator*(const Matrix& o) const {
        if (cols != o.rows) throw runtime_error("Columns of A must equal rows of B for multiplication.");
        Matrix r(rows, o.cols);
        for (int i = 0; i < rows; i++)
            for (int k = 0; k < cols; k++)
                for (int j = 0; j < o.cols; j++) r.a[i][j] += a[i][k] * o.a[k][j];
        return r;
    }

    Matrix transpose() const {
        Matrix t(cols, rows);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) t.a[j][i] = a[i][j];
        return t;
    }

    // Determinant via Gaussian elimination with partial pivoting: O(n^3)
    double determinant() const {
        if (rows != cols) throw runtime_error("Determinant requires a square matrix.");
        int n = rows;
        Matrix m = *this;
        double det = 1.0;
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int i = col + 1; i < n; i++)
                if (fabs(m.a[i][col]) > fabs(m.a[pivot][col])) pivot = i;
            if (fabs(m.a[pivot][col]) < EPS) return 0.0;
            if (pivot != col) {
                swap(m.a[pivot], m.a[col]);
                det = -det;  // a row swap flips the sign
            }
            det *= m.a[col][col];
            for (int i = col + 1; i < n; i++) {
                double f = m.a[i][col] / m.a[col][col];
                for (int j = col; j < n; j++) m.a[i][j] -= f * m.a[col][j];
            }
        }
        return det;
    }

    // Gauss-Jordan on [A | I]. Returns nullopt if singular.
    optional<Matrix> inverse() const {
        if (rows != cols) throw runtime_error("Inverse requires a square matrix.");
        int n = rows;
        Matrix m = *this, inv = identity(n);
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int i = col + 1; i < n; i++)
                if (fabs(m.a[i][col]) > fabs(m.a[pivot][col])) pivot = i;
            if (fabs(m.a[pivot][col]) < EPS) return nullopt;
            swap(m.a[pivot], m.a[col]);
            swap(inv.a[pivot], inv.a[col]);
            double p = m.a[col][col];
            for (int j = 0; j < n; j++) {
                m.a[col][j] /= p;
                inv.a[col][j] /= p;
            }
            for (int i = 0; i < n; i++) {
                if (i == col) continue;
                double f = m.a[i][col];
                for (int j = 0; j < n; j++) {
                    m.a[i][j] -= f * m.a[col][j];
                    inv.a[i][j] -= f * inv.a[col][j];
                }
            }
        }
        return inv;
    }

    // Reduced row echelon form; also reports the rank (number of pivots).
    Matrix rref(int* rankOut = nullptr) const {
        Matrix m = *this;
        int lead = 0;  // current pivot row
        for (int col = 0; col < cols && lead < rows; col++) {
            int pivot = lead;
            for (int i = lead + 1; i < rows; i++)
                if (fabs(m.a[i][col]) > fabs(m.a[pivot][col])) pivot = i;
            if (fabs(m.a[pivot][col]) < EPS) continue;  // no pivot in this column
            swap(m.a[pivot], m.a[lead]);
            double p = m.a[lead][col];
            for (int j = 0; j < cols; j++) m.a[lead][j] /= p;
            for (int i = 0; i < rows; i++) {
                if (i == lead) continue;
                double f = m.a[i][col];
                for (int j = 0; j < cols; j++) m.a[i][j] -= f * m.a[lead][j];
            }
            lead++;
        }
        if (rankOut) *rankOut = lead;
        return m;
    }

    int rank() const {
        int r = 0;
        rref(&r);
        return r;
    }
};

void menu() {
    cout << "\n===== Matrix Calculator =====\n"
         << " 1. Add            (A + B)\n"
         << " 2. Subtract       (A - B)\n"
         << " 3. Multiply       (A * B)\n"
         << " 4. Transpose      (A)\n"
         << " 5. Determinant    (A)\n"
         << " 6. Inverse        (A)\n"
         << " 7. Row-reduce     (RREF of A)\n"
         << " 8. Rank           (A)\n"
         << " 0. Exit\n"
         << "Choice: ";
}

int main() {
    int choice;
    while (true) {
        menu();
        if (!(cin >> choice)) break;
        if (choice == 0) break;
        if (choice < 1 || choice > 8) {
            cout << "Invalid choice.\n";
            continue;
        }
        try {
            Matrix A, B;
            A.read("A");
            if (choice <= 3) B.read("B");

            switch (choice) {
                case 1: (A + B).print("A + B ="); break;
                case 2: (A - B).print("A - B ="); break;
                case 3: (A * B).print("A * B ="); break;
                case 4: A.transpose().print("A^T ="); break;
                case 5: cout << "det(A) = " << fixed << setprecision(4) << A.determinant() << "\n"; break;
                case 6: {
                    auto inv = A.inverse();
                    if (inv) inv->print("A^-1 =");
                    else cout << "Matrix is singular (determinant = 0); no inverse exists.\n";
                    break;
                }
                case 7: A.rref().print("RREF(A) ="); break;
                case 8: cout << "rank(A) = " << A.rank() << "\n"; break;
            }
        } catch (const exception& e) {
            cout << "Error: " << e.what() << "\n";
        }
    }
    cout << "Goodbye!\n";
    return 0;
}
