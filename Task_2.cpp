#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
vector<double> Gauss(vector<vector<double>> matrix_A, vector<double> vector_b, int n)
{
	int p;
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

int main()
{
	setlocale(LC_ALL, "Russian");
	int n;
	cout << "Задание 2. Итерационные методы решения линейных систем.\n"
		<< "\nВведите размерность исходной матрицы:\n"
		<< "n = ";
	cin >> n;
	vector<vector<double>> matrix_A(n, vector<double>(n));
	vector<double>vector_b(n), exSolution(n); //exact solution

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
	exSolution = Gauss(matrix_A, vector_b, n);
	cout << "\nРешение методом Гаусса:\n";
	for (int i = 0; i < n; i++)
		cout << "x[" << i + 1 << "] = " << exSolution[i] << "\n";

	return(0);
}