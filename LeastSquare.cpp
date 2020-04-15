#include <iostream> 
#include <iomanip>
#include <cmath>	
#include <string>
#include <cstdio>

#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

using namespace std;
class Matrix {
public:
	int n, m;
	double** storage;
	Matrix(int a, int b) {
		n = a;
		m = b;
		storage = new double* [a];
		for (int i = 0; i < n; i++) {
			storage[i] = new double[b];
		}
	}
	friend ostream& operator << (ostream& out, const Matrix& mat);
	friend istream& operator >> (istream& in, Matrix& mat);
	Matrix operator = (Matrix const& other) {
		Matrix rez(other.n, other.m);
		n = other.n, m = other.m;
		for (int i = 0; i < other.n; i++) {
			for (int j = 0; j < other.m; j++) {
				rez.storage[i][j] = other.storage[i][j];
				storage[i][j] = other.storage[i][j];
			}
		}
		return rez;
	}
	Matrix operator + (Matrix const& other) {
		Matrix rez(other.n, other.m);
		for (int i = 0; i < other.n; i++) {
			for (int j = 0; j < other.m; j++) {
				rez.storage[i][j] = storage[i][j] + other.storage[i][j];
			}
		}
		return rez;
	}
	Matrix operator - (Matrix const& other) {
		Matrix rez(other.n, other.m);
		for (int i = 0; i < other.n; i++) {
			for (int j = 0; j < other.m; j++) {
				rez.storage[i][j] = storage[i][j] - other.storage[i][j];
			}
		}
		return rez;
	}
	Matrix operator * (Matrix const& other) {
		Matrix rez(n, other.m);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < other.m; j++) {
				rez.storage[i][j] = 0;
				for (int k = 0; k < m; k++) {
					rez.storage[i][j] += storage[i][k] * other.storage[k][j];
				}
			}
		}
		return rez;
	}
	Matrix T() {
		Matrix rez(m, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				rez.storage[j][i] = storage[i][j];
			}
		}
		return rez;
	}
};


class squareMatrix : public Matrix {
public:
	squareMatrix(int n) : Matrix(n, n) {};
};

class identityMatrix : public squareMatrix {
public:
	identityMatrix(int n) :squareMatrix(n) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					storage[i][j] = 1.00;
				}
				else {
					storage[i][j] = 0.00;
				}
			}
		}
	};

};

class eliminationMatrix : public identityMatrix {
public:
	eliminationMatrix(int n, int i, int j, double val) : identityMatrix(n) {
		storage[i][j] = val;
	}
};

class permutationMatrix : public identityMatrix {
public:
	permutationMatrix(int n, int i, int j) :identityMatrix(n) {
		double* temp = new double[n];
		for (int k = 0; k < n; k++) {
			temp[k] = storage[j][k];
		}
		for (int k = 0; k < n; k++) {
			storage[j][k] = storage[i][k];
		}
		for (int k = 0; k < n; k++) {
			storage[i][k] = temp[k];
		}
	}
};

ostream& operator << (ostream& out, const Matrix& mat)
{
	double e = 0.000000001;
	for (int i = 0; i < mat.n; i++) {
		for (int j = 0; j < mat.m; j++) {
			if (j == mat.m - 1) { out << setprecision(2) << fixed << mat.storage[i][j]+e; }
			else { out << setprecision(2) << fixed << mat.storage[i][j]+e << " "; }
		}
		out << endl;
	}
	return out;
}

istream& operator >> (istream& in, Matrix& mat)
{
	for (int i = 0; i < mat.n; i++) {
		for (int j = 0; j < mat.m; j++) {
			in >> mat.storage[i][j];
		}
	}
	return in;
}

class CV : public Matrix {
public:
	CV(int n) : Matrix(n, 1) {};
};

Matrix inverse(int n, Matrix p) {
	squareMatrix a(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a.storage[i][j] = p.storage[i][j];
		}
	}

	identityMatrix left(n);
	for (int j = 0; j < n - 1; j++) {
		for (int i = j + 1; i < n; i++) {
			double val = a.storage[i][j] / a.storage[j][j];
			eliminationMatrix e(n, i, j, -val);
			(Matrix)a = e * a;
			(Matrix)left = e * left;
		}
	}
	for (int j = n - 1; j > 0; j--) {
		for (int i = j - 1; i > -1; i--) {
			if (a.storage[i][j] != 0) {
				double val = a.storage[i][j] / a.storage[j][j];
				eliminationMatrix e(n, i, j, -val);
				(Matrix)a = e * a;
				(Matrix)left = e * left;
			}
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			left.storage[i][j] /= a.storage[i][i];
		}
		a.storage[i][i] = 1;
	}
	return left;

}

int main() {
	int m;
	cin >> m;
	CV b(m);
	double* t = new double[m];
	for (int i = 0; i < m; i++) {
		cin >> t[i] >> b.storage[i][0];
	}
	int n;
	cin >> n;
	Matrix A(m,n+1);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n+1; j++) {
			A.storage[i][j] = pow(t[i], j);
		}
	}
	cout << "A:" << endl;
	cout << A;
	cout << "A_T*A:" << endl;
	Matrix T = A.T() * A;
	cout << T;
	cout << "(A_T*A)^-1:" << endl;
	Matrix I = inverse(n+1, T);
	cout << I;
	Matrix B = A.T() * b;
	cout << "A_T*b:" << endl << B;
	cout << "x~:" << endl << I * B;
	Matrix X = I * B;

	

	FILE* pipe = _popen(GNUPLOT_NAME, "w");
	if (pipe != NULL) {
		string s = "f(x)=";
		for (int j = 0; j < n + 1; j++) {
			s.append("(");
			s.append(to_string(X.storage[j][0]));
			s.append("*x**");
			s.append(to_string(j));
			s.append(")");
			if (j != n) {
				s.append("+");
			}
		}
		s.append("\n");
		fprintf(pipe, "set term wx\n");
		fprintf(pipe, s.c_str());
		fprintf(pipe, "plot '-' using 1:2 title 'exp' with points pointtype 5 pointsize 1,  f(x) title 'appr' with lines\n");
		for (int i = 0; i < m; i++) {
			fprintf(pipe, "% f % f\n", t[i], b.storage[i][0]);

		}
		fprintf(pipe, "%s\n", "e");
		fflush(pipe);


		_pclose(pipe);
		

	}
	else {
		std::cout << "Could not open pipe" << std::endl;
	}

	return 0;
}