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
double normVector(vector<double> vector)
{
	double norm = 0;
	for (int i = 0; i < vector.size(); i++)
		norm += fabs(vector[i]);
	return norm;
}
double normMatrix(vector<vector<double>>matrix)
{
	int n = matrix.size();
	vector<double>sums(n);
	for (int j = 0; j < n; j++)
	{
		double sum = 0;
		for (int i = 0; i < n; i++)
			sum += fabs(matrix[i][j]);
		sums[j] = sum;
	}
	double max = *max_element(sums.begin(), sums.end());
	return max;
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
			if (abs(a[p][k]) < abs(a[m][k])) //поиск максимального ведущего элемента
				p = m;
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
				a[m][i] = a[m][i] - c * a[k][i];
		}
	}
	x[n - 1] = b[n - 1] / a[n - 1][n - 1];
	for (int k = n - 1; k >= 0; k--)
	{
		s = 0;
		for (int i = k + 1; i < n; i++)				//обратный ход метода Гаусса
			s = s + a[k][i] * x[i];
		x[k] = (b[k] - s) / a[k][k];
	}
	return x;
}
vector<vector<double>> matrix_D(vector<vector<double>>matrix_A)
{
	int n = matrix_A.size();
	vector<vector<double >> D(n,vector<double>(n));
	for(int i=0; i<n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				D[i][i] = matrix_A[i][i];
			else
				D[i][j] = 0;
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
	//show(multiply);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				H_D[i][j] = 1 - multiply[i][j];
			else
				H_D[i][j] = 0 - multiply[i][j];
		}
	return H_D;
}
int priori_estimate_step(double norm_H, vector<double>g)
{
	int n = g.size();
	int k = 0;
	double sum = 0, estimate = 1;
	double norm_g = normVector(g);
	//estimate = ((pow(norm_H, p) / (1 - norm_H)) * norm_g);
	while (estimate > 0.001)
	{
		estimate = ((pow(norm_H, k) / (1 - norm_H)) * norm_g);
		k++;
	}
	//return estimate;
	return k;
}
vector<vector<double>> simple_itteration(vector<vector<double>> alpha, vector<double>beta, vector<double>exSol)
{
	int n = beta.size(),step = 0;
	double e = 0.001, factEstimate=1;
	vector<double> x_k(n), x_0(n),estimate(n),Ax(n);
	vector<vector<double>> matrix(3, vector<double>(n));
	for (int i = 0; i < n; i++)
		x_0[i] = beta[i];
	while (factEstimate > e)
	{
		Ax = multiplyMatrixVector(alpha, x_0);
		for (int i = 0; i < n; i++)
		{
			x_k[i] = Ax[i] + beta[i];
			estimate[i] = x_k[i] - exSol[i];
			for (int m = 0; m < n; m++)
			{
				matrix[0][m] = x_k[m];
				matrix[1][m] = x_0[m];
			}
			x_0[i] = x_k[i];
		}
		factEstimate = normVector(estimate);
		step++;
	}
	matrix[2][0] = step;
	return matrix;
}
double apriori_estimate(vector<vector<double>> H, vector<double> g, int k)
{
	double normM = normMatrix(H), normV = normVector(g);
	return normM * normV + (pow(normM, k) / (1 - normM)) * normV;
}
double aposteriori_estimate(vector<vector<double>> H, vector<vector<double>> x)
{
	double normM = normMatrix(H);
	vector<double> delta(H.size());
	for (int i = 0; i < H.size(); i++)
		delta[i] = x[0][i] - x[1][i];
	double normD = normVector(delta);
	return (normM/(1-normM))*normD;
}
double scalar(vector<double> a, vector<double> b)
{
	double ab = 0;
	for (int i = 0; i < a.size(); i++)
		ab += a[i] * b[i];
	return ab;
}
double spectr_radius(vector<vector<double>> matrix)
{
	int n = matrix.size();
	double lambda = 0;
	vector<double> x_0(n), x_k(n);
	for (int i = 0; i < n; i++)
		x_0[i] = 1;
	for (int k = 0; k < 50; k++)
	{
		x_k = multiplyMatrixVector(matrix, x_0);
		lambda = scalar(x_k, x_0) / scalar(x_0, x_0);
		x_0 = x_k;
	}
	return fabs(lambda);
}
vector<double> Lusterink(vector<vector<double>> H, vector<vector<double>> solve)
{
	int n = H.size();
	vector<double> x_lust(n);
	for (int i = 0; i < n; i++)
		x_lust[i] = solve[1][i] + (1 / (1 - spectr_radius(H))) * (solve[0][i] - solve[1][i]);
	return x_lust;
}
vector<double> Zeydel(vector<vector<double>>H, vector<double> g, vector<double> exSol)
{
	int n = g.size(), step = 0;
	double e = 0.001, R = 1, sum1=0, sum2=0;
	vector<double>x_k(n), x_0(n),delta(n);
	for (int i = 0; i < n; i++)
	{
		x_0[i] = g[i];
		x_k[i] = 0;
	}

	while (R > e)
	{
		for (int i = 0; i < n;i++)
		{
			sum1 = 0; sum2 = 0;
			for (int j = 0; j <= i - 1; j++)
				sum1 += H[i][j] * x_k[j];
			for (int j = i; j < n; j++)
				sum2 += H[i][j] * x_0[j];
			x_k[i] = sum1 + sum2 + g[i];
		}
		for (int i = 0; i < n; i++)
			delta[i] = fabs(x_k[i] - exSol[i]);
		R = normVector(delta);
		x_0 = x_k;
		step++;
	}
	x_k.push_back(step);
	return x_k;
}
vector<double> relax(vector<vector<double>> H, vector<double> g, vector<double> exSol)
{
	int n = g.size(), step = 0;
	double R = 1, e = 0.001, sum = 0, sum1, sum2;
	vector<double> x_0(n), x_k(n+1),delta(n);
	double q = 2 / (1 + sqrt(1 - pow(spectr_radius(H), 2)));
	for (int i = 0; i < n; i++)
	{
		x_0[i] = 1;
		x_k[i] = 0;
	}	
	while (R >= e)
	{
		for (int i = 0; i < n; i++)
		{
			sum = 0; sum1 = 0; sum2 = 0;
			for (int j = 0; j <= i - 1; j++)
				sum1 += H[i][j] * x_k[j];
			for (int j = i+1; j < n; j++)
				sum2 += H[i][j] * x_0[j];
			sum = sum1 + sum2 - x_0[i] + g[i];
			x_k[i] = x_0[i] + q * sum;
		}
		x_0 = x_k;
		for (int i = 0; i < n; i++)
			delta[i] = (exSol[i] - x_k[i]);
		R = normVector(delta);		
		step++;
	}
	x_k[n] = step;
	return x_k;
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/

int main()
{
	setlocale(LC_ALL, "Russian");
	int n;
	double norm_HD = 0;
	cout << "Задание 2. Итерационные методы решения линейных систем. Вариант 3.\n"
		<< "\nВведите размерность исходной матрицы:\n"
		<< "n = ";
	cin >> n;
	vector<vector<double>> matrix_A(n, vector<double>(n)),
		D(n, vector<double>(n)),
		reverse_D(n, vector<double>(n)),
		H_D(n, vector<double>(n)), solve_simple_itteration(3, vector<double> (n));

	vector<double>vector_b(n),
		exSolution(n),//exact solution
		g_D(n),
		estimate(n),
		x_lust(n),
		x_zeydel(n+1),
		x_relax(n+1),
		delta(n);

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
	//(1)
	exSolution = Gauss(matrix_A, vector_b);
	cout << "\nРешение методом Гаусса:\n";
	for (int i = 0; i < n; i++)
		cout << "x[" << i + 1 << "] = " << exSolution[i] << "\n";
	//(2)
	cout << "\nМатрица D:\n";
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
	cout << "\n\nМатрица H_D:\n";
	show(H_D);
	norm_HD = normMatrix(H_D);
	cout << "||H_D|| = " << norm_HD << "\n";
	//(3)
	cout << "Априорная оценка при k = " << priori_estimate_step(norm_HD, g_D)<<"\n";

	//проверка на сходимость
	if (norm_HD > 1)
	{
		cout << "\nУсловие на сходимость не выполняется.\n";
		return(0);
	} 
	//(4)
	solve_simple_itteration = simple_itteration(H_D, g_D, exSolution);
	cout << "\nРешение методом простой итерации:\n";
	for (int i = 0; i < n; i++)
		cout <<"x["<<i+1<<"] = " << solve_simple_itteration[0][i] << "\n";
	int step = solve_simple_itteration[2][0];
	cout << "Фактическое число итераций k = " << step << "\n";
	cout << "Априорная оценка = " << apriori_estimate(H_D, g_D, step)<<"\n";
	cout << "Апостериорная оценка = " << aposteriori_estimate(H_D, solve_simple_itteration)<<"\n";
	x_lust = Lusterink(H_D, solve_simple_itteration);
	cout << "Последнее приближение по Люстеринку:\n";
	for (int i = 0; i < n; i++)
		cout << "x[" << i + 1 << "] = " << x_lust[i] << "\n";
	for (int i = 0; i < n; i++)
		delta[i] = x_lust[i] - exSolution[i];
	cout << "Фактическая погрешность приближения по Люстеринку = " << normVector(delta);
	//(5)
	x_zeydel = Zeydel(H_D, g_D, exSolution);
	cout << "\n\nРешение методом Зейделя:\n";
	for (int i = 0; i < n; i++)
		cout << "x[" << i + 1 << "] = " << x_zeydel[i] << "\n";
	step = x_zeydel[n];
	cout << "Фактическое число итераций k = " << step << "\n";
	//(6)
	cout << "\nСпектральный радиус матрицы перехода = " << spectr_radius(H_D)<<"\n";
	//(7)
	cout << "\nРешение методом верхней релаксации:\n";
	x_relax = relax(H_D, g_D, exSolution);
	for (int i = 0; i < n; i++)
		cout << "x[" << i + 1 << "] = " << x_relax[i] << "\n";
	step = x_relax[n];
	cout << "Фактическое число итераций k = " << step << "\n";
	return(0);
}