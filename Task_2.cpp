#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

void show(vector<vector<double>>matrix)
{
	int n = matrix.size();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << "\t" << matrix[i][j];
		cout << endl;
	}
	cout << endl;
}
vector<double> multiplyMatrixVector(vector<vector<double>>A, vector<double>b)
{
	int n = A.size();
	vector<double>Ab(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Ab[i] += A[i][j] * b[j];
	return Ab;
}
vector<vector<double>> multiplyMatrixMatrix(vector <vector <double>> A, vector <vector <double>> B)
{
	int n = A.size();
	vector <vector <double>>AB(n,vector<double>(n));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				AB[i][j] += A[i][k] * B[k][j];
	return AB;
}
vector<double> Gauss(vector<vector<double>> matrix_A, vector<double> vector_b)
{
	int p, n = matrix_A.size();
	double r, c, s;
	vector<double>x(n), b(n);
	vector<vector<double>>a(n, vector<double>(n));
	a = matrix_A;
	b = vector_b;
	for (int k = 0; k < n; k++)
	{
		p = k;
		for (int m = k + 1; m < n; m++)
		{
			if (abs(a[p][k]) < abs(a[m][k])) //поиск максимального ведущего элемента
			{
				p = m;
			}
		}
		for (int j = k; j < n; j++)
		{
			r = a[k][j];
			a[k][j] = a[p][j];   //перестановка строк
			a[p][j] = r;
		}
		r = b[k];
		b[k] = b[p];   //перестановка свободных членов
		b[p] = r;
		for (int m = k + 1; m < n; m++)
		{
			c = a[m][k] / a[k][k];
			b[m] = b[m] - c * b[k]; //приведение матрицы к верхнетреугольному виду
			for (int i = k; i < n; i++)
			{
				a[m][i] = a[m][i] - c * a[k][i];
			}
		}
	}
	x[n - 1] = b[n - 1] / a[n - 1][n - 1];
	for (int k = n - 1; k >= 0; k--)
	{
		s = 0;
		for (int i = k + 1; i < n; i++)				//обратный ход метода Гаусса
		{
			s = s + a[k][i] * x[i];
		}
		x[k] = (b[k] - s) / a[k][k];
	}
	return x;
}
vector<vector<double>> matrix_D(vector<vector<double>>matrix_A)
{
	int n = matrix_A.size();
	vector<vector<double >> D(n,vector<double>(n));
	for(int i=0; i<n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				D[i][i] = matrix_A[i][i];
			else
				D[i][j] = 0;
		}
	}
	return D;
}
vector<vector<double>> reverse_matrix(vector<vector<double>> matrix)
{
	int n = matrix.size();
	vector<vector<double>> reverse(n, vector<double>(n));
	vector<double>b(n);
	vector<double>solve(n);
	for (int i = 0; i < n; i++)
	{
		b[i] = 1;
		solve = Gauss(matrix, b);
		for (int j = 0; j < n; j++)
			reverse[j][i] = solve[j];
		b[i] = 0;
	}
	return reverse;
}
vector<vector<double>> matrix_HD(vector<vector<double>>matrix, vector<vector<double>>reverse)
{
	int n = matrix.size();
	vector<vector<double>>multiply(n, vector<double>(n)), 
		H_D(n, vector<double > (n));
	multiply = multiplyMatrixMatrix(reverse, matrix);
	show(multiply);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			H_D[i][j] = 1 - multiply[i][j];
	return H_D;
}
double normMatrix(vector<vector<double>>matrix)
{
	int n = matrix.size();
	vector<double>sums(n);
	for (int j = 0; j < n; j++)
	{
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += fabs(matrix[i][j]);
		}
		sums[j] = sum;
	}
	double max = *max_element(sums.begin(), sums.end());
	return max;
}
int priori_estimate(double norm_H,vector<double>&g )
{
	int n = g.size(),k=0;
	double norm_g = 0,sum=0, estimate = 1;
	for (int i = 0; i < n; i++)
		sum += g[i] * g[i];
	norm_g = sqrt(sum);
	while (estimate > 0.001)
	{
		estimate = ((pow(norm_H, k) / (1 - norm_H))* norm_g);
		k++;
	}
	return k;
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/

int main()
{
	setlocale(LC_ALL, "Russian");
	int n;
	double e = 0.001, norm_HD=0;
	cout << "Задание 2. Итерационные методы решения линейных систем. Вариант 3.\n"
		<< "\nВведите размерность исходной матрицы:\n"
		<< "n = ";
	cin >> n;
	vector<vector<double>> matrix_A(n, vector<double>(n)), 
		D(n,vector<double>(n)),
		reverse_D(n,vector<double>(n)),
		H_D(n,vector<double>(n));

	vector<double>vector_b(n), 
		exSolution(n),
		g_D(n); //exact solution

	cout << "\nЗаполняем матрицу (А):\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			{
				cout << "[" << i + 1 << "][" << j + 1 << "] = ";
				cin >> matrix_A[i][j];
			}
		}
	}

	cout << "\nЗаполняем вектор (b):\n";
	for (int i = 0; i < n; i++)
	{
		cout << "b[" << i + 1 << "] = ";
		cin >> vector_b[i];
	}

	cout << "\nРасширенная матрица (А|b):\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << "\t" << matrix_A[i][j];
		cout << "\t" << vector_b[i] << endl;
	}

	exSolution = Gauss(matrix_A, vector_b);
	cout << "\nРешение методом Гаусса:\n";
	for (int i = 0; i < n; i++)
		cout << "x[" << i + 1 << "] = " << exSolution[i] << "\n";

	cout << "Матрица D:\n";
	D = matrix_D(matrix_A);
	show(D);

	cout << "Обратная матрица D^(-1):\n";
	reverse_D = reverse_matrix(D);
	show(reverse_D);

	g_D = multiplyMatrixVector(reverse_D, vector_b);
	H_D = matrix_HD(matrix_A, reverse_D);
	cout << "Вектор g_D:\n";
	for (int i = 0; i < n; i++)
		cout << "\t" << g_D[i];
	cout << "\nМатрица H_D:\n";
	show(H_D);
	norm_HD = normMatrix(H_D);
	cout << "||H_D|| = " << norm_HD << "\n";
	cout << "Априорная оценка k = " << priori_estimate(norm_HD, g_D);

		
	
	return(0);
}