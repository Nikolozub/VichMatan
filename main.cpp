#include <iostream>
#include <conio.h>
#include <algorithm>

using namespace std;

void swap_colums( int size, double **matrix, int c_1, int c_2 )
{
	if ( ( size <= 0 ) || ( c_1 < 0) || ( c_1 >= size)  || ( c_2 < 0) || ( c_2 >= size) ) return;

	for ( int i = 0; i < size; i++ )
		swap( matrix[i][c_1], matrix[i][c_2] );
}

void printMatrix( int n, int m, double **matrix )
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

// Метод Гаусса с выбором главного элемента по строке 
// Возврат результата через параметр x_matrix
int gauss( int size, double **a_matrix, double *b_matrix, double *x_matrix )
{
	if ( size <= 0 ) return 0;

	printMatrix(size, size, a_matrix);
	cout << endl;
	printMatrix( 1, size, &b_matrix );
	cout << endl;

	for ( int k = 0; k < size - 1; k++ )
	{
		cout << "k = " << k << endl;

		int max = k;
		for ( int i = k + 1; i < size; i++ ) 
		{
			if ( abs( a_matrix[k][i] ) > abs( a_matrix[k][max] ) ) max = i;
		}

		if (max != k) swap_colums( size, a_matrix, k, max );
	    printMatrix( size, size, a_matrix );
	    cout << endl;

		for ( int m = k + 1; m < size; m++ )
		{
			double koef = a_matrix[m][k] / a_matrix[k][k];

			b_matrix[m] = b_matrix[m] - koef * b_matrix[k];

			for ( int l = k; l < size; l++ )
			{
				a_matrix[m][l] = a_matrix[m][l] - koef * a_matrix[k][l];
			}			
		}

	    printMatrix( size, size, a_matrix );
	    cout << endl;
		printMatrix( 1, size, &b_matrix );
	    cout << endl;
	}


	for ( int i = size - 1; i >= 0; i-- )
	{
		double sum = 0;
		for ( int j = i + 1; j < size; j++ )
		{
			sum += a_matrix[i][j] * x_matrix[j];
		}

		x_matrix[i] = ( b_matrix[i] - sum ) / a_matrix[i][i];
	}

	return 1;
}

int norma1(int size, double **matrix)
{	
	int max = 0;

	for (int i = 0; i < size; i++)
	{
		int sum = 0;

		for (int j = 0; j < size; j++)
		{
			sum += abs(matrix[i][j]);
		}

		if (sum > max) max = sum;
	}

	return max;
}

bool shodimost(int size, double **matrix)
{
	int sum = 0;

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < i; j++)
			sum += abs(matrix[i][j]);
		for (int j = i + 1; j < size; j++)
			sum += abs(matrix[i][j]);

		if (matrix[i][i] <= sum) return false;
		sum = 0;
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < i; j++)
			sum += abs(matrix[j][i]);
		for (int j = i + 1; j < size; j++)
			sum += abs(matrix[j][i]);

		if (matrix[i][i] <= sum) return false;
		sum = 0;
	}

	return true;
}

// Метод Зейделя
// Возврат результата через параметр x_matrix
int zeidel ( int size, double **a_matrix, double *b_matrix, double e, double *x_matrix )
{
	if ( size <= 0 ) return 0;

	for ( int i = 0; i < size; i++ )
		x_matrix[i] = 0;
	
	for (int q = 0; q < 3; q++)
	{
		for ( int i = 0; i < size; i++ )
		{
			x_matrix[i] = b_matrix[i];

			for ( int j = 0; j < size; j++ )
			{
				if (i != j) x_matrix[i] -= a_matrix[i][j] * x_matrix[j];
			}
			x_matrix[i] /= a_matrix[i][i];
		}
	}

	return 1;
}


int main()
{
	const int n = 3;

	double A[n][n] = {
		1 , 9 ,  1 ,
		2 , 2 , 11 ,
	    10 , 2 , 1
	};

	double B[n] = { 12, 16, 14 };

	double *X = new double[n];

	double **A2 = new double*[n];
	for (int i = 0; i < n; i++)
		A2[i] = new double[n];

	A2[0][0] = 1;
	A2[0][1] = 9;
	A2[0][2] = 1;

	A2[1][0] = 2;
	A2[1][1] = 2;
	A2[1][2] = 11;

	A2[2][0] = 10;
	A2[2][1] = 2;
	A2[2][2] = 1;

	//

	//gauss(n, A2, B, X);
	zeidel(n, A2, B, 10, X);
	printMatrix( 1, n, &X );

	cout << endl << shodimost(n, A2);

	_getch();
	return 0;
}