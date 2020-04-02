#include <iostream>
#include <conio.h>
#include <algorithm>
#include <limits>

#define EPSILON 0.00000001

using namespace std;

typedef std::numeric_limits< double > dbl;

// меняет местами столбцы квадратной матрицы
void swap_colums( int size, double **matrix, int c_1, int c_2 )
{
	if ( ( size <= 0 ) || ( c_1 < 0) || ( c_1 >= size)  || ( c_2 < 0) || ( c_2 >= size ) || ( c_1 == c_2 ) ) return;

	for ( int i = 0; i < size; i++ )
		swap( matrix[i][c_1], matrix[i][c_2] );
}

// печатает матрицу на экран
void printMatrix( int n, int m, const double * const *matrix )
{
	for ( int i = 0; i < n; i++ )
	{
		for ( int j = 0; j < m; j++ )
		{
			cout << matrix[i][j] << '\t';
		}

		cout << '\n';
	}
}

// копирует матрицу src в dest
void copyMatrix( int n, int m, const double * const *src, double **dest )
{
	for ( int i = 0; i < n; i++ ) 
		for ( int j = 0; j < m; j++ ) 
				dest[i][j] = src[i][j];		
}


// Метод Гаусса с выбором главного элемента по строке 
// Возврат результата через параметр X_matrix
int gauss( int size, const double * const *A_matrix, const double *B_matrix, double *X_matrix )
{
	if ( size <= 0 ) return 0;

	// копируем матрицы
	double **a_matrix = new double*[size];
	for (int i = 0; i < size; i++)
			a_matrix[i] = new double[size];
	copyMatrix( size, size, A_matrix, a_matrix );

	double *b_matrix = new double[size];
	copyMatrix(1, size, &B_matrix, &b_matrix);

	// порядок следования X
	int *order = new int[size];
	for (int i = 0; i < size; i++)
			order[i] = i;

	// формируем треугольную матрицу с перестановкой столбцов
	for ( int k = 0; k < size; k++ )
	{

		int max = k;
		for ( int i = k + 1; i < size; i++ ) 
		{
			if ( abs( a_matrix[k][i] ) > abs( a_matrix[k][max] ) ) max = i;
		}

		if ( a_matrix[k][max] < EPSILON ) 
		{
			for ( int v = 0; v < size; v++ ) 
				delete[] a_matrix[v];
			delete[] a_matrix;
			delete[] b_matrix;
			delete[] order;
			return 0;
		}

		if (max != k) 
		{
			swap_colums( size, a_matrix, k, max );
			swap(order[k], order[max]);
		}


		for ( int m = k + 1; m < size; m++ )
		{
			double koef = a_matrix[m][k] / a_matrix[k][k];

			b_matrix[m] = b_matrix[m] - koef * b_matrix[k];

			a_matrix[m][k] = 0.0;

			for ( int l = k + 1; l < size; l++ )
			{
				a_matrix[m][l] = a_matrix[m][l] - koef * a_matrix[k][l];
			}			
		}

	}

	// расчитываем X с учетом их порядка следования
	for ( int i = size - 1; i >= 0; i-- )
	{
		double sum = 0;
		for ( int j = i + 1; j < size; j++ )
		{
			sum += a_matrix[i][j] * X_matrix[ order[j] ];
		}

		X_matrix[ order[i] ] = ( b_matrix[i] - sum ) / a_matrix[i][i];
	}

	for ( int v = 0; v < size; v++ ) 
		delete[] a_matrix[v];
	delete[] a_matrix;
	delete[] b_matrix;
	delete[] order;

	return 1;
}


// Метод Зейделя
// Возврат результата через параметр x_matrix
int zeidel ( int size, const double * const *a_matrix, const double *b_matrix, double e, double *x_matrix )
{
	if ( size <= 0 ) return 0;

	double sum = 0;
	double mNorma = 0;

	int i;

	for (i = 0; i < size; i++)
	{
		for (int j = 0; j < i; j++)
			sum += abs(a_matrix[i][j]);
		for (int j = i + 1; j < size; j++)
			sum += abs(a_matrix[i][j]);

		sum /= abs( a_matrix[i][i] );

		if ( sum >= 1 ) return 0;

		if ( sum > mNorma ) mNorma = sum;

		sum = 0;
	}

	// начальное приближение
	for ( int i = 0; i < size; i++ )
		x_matrix[i] = 0;

	double limit = ( ( 1 - mNorma ) / mNorma ) * e;

	double dx_norma;

	double prevX, curX;


	do
	{
		dx_norma = 0;

		for ( int i = 0; i < size; i++ )
		{
			prevX = x_matrix[i];

			curX = b_matrix[i];

			for ( int j = 0; j < size; j++ )
			{
				if (i != j) curX -= a_matrix[i][j] * x_matrix[j];
			}
			curX /= a_matrix[i][i];

			dx_norma += abs( curX - prevX );
			x_matrix[i] = curX;

		}

	} while ( dx_norma >= limit );

	return 1;
}


int main()
{
	const int n = 3;
	const double e = 0.0001;

	double **A = new double*[n];
	for (int i = 0; i < n; i++)
		A[i] = new double[n];

	A[0][0] = 10;
	A[0][1] = 1;
	A[0][2] = 2;

	A[1][0] = 2;
	A[1][1] = 11;
	A[1][2] = 2;

	A[2][0] = 1;
	A[2][1] = 1;
	A[2][2] = 9;

	double *B = new double[n];
	
	B[0] = 14;
	B[1] = 26;
	B[2] = 12;

	double *X = new double[n];

	cout.precision( dbl::max_digits10 );

	cout << "A[][]:\n";
	printMatrix( n, n, A );
	cout << "\nB[]:\n";
	printMatrix( 1, n, &B );


	cout << "\n\nGauss method: \n";
	if ( gauss( n, A, B, X ) ) 
	{
		cout << "X[]: ";
		printMatrix( 1, n, &X );
	}
	else
	{
		cout << "error\n";
	}

	cout << "\nZeidel method: ( e = " << e << " )\n";
	if ( zeidel( n, A, B, e, X ) ) 
	{
		cout << "X[]: ";
		printMatrix( 1, n, &X );
	}
	else
	{
		cout << "error\n";
	}

	_getch();
	return 0;
}