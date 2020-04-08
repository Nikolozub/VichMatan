using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace priplijenie
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        // интерполяция многочленом в форме лагранжа
        private double lagrange( double[] X, double[] Y, double x ) 
        {
            double y = 0;

            double k;
            for (int i = 0; i < X.Length; i++)
            {
                k = Y[i];
                for (int j = 0; j < X.Length; j++)
                {
                    if (i != j) k *= (x - X[j]) / (X[i] - X[j]);    
                }
                y += k;
            }

            return y;
        }

        double razd_raznost(double[] X, double[] Y, int i, int n) 
        {
            if (n == i) return Y[i];
            return (razd_raznost(X, Y, i, n-1) - razd_raznost(X, Y, i+1, n)) / (X[i] - X[n]);
        }

        double newton(double[] X, double[] Y, int k, double x) 
        {
            if (k >= X.Length) return 0;
            return razd_raznost(X, Y, 0, k) + (x - X[k]) * newton(X, Y, k + 1, x);
        }

        double calc_polynom(double[] A, double x) 
        {
            double y = 0;
            for (int i = A.Length - 1; i >= 0; i--)
            {
                y = A[i] + y * x;
            }
            return y;
        }

        // меняет местами столбцы квадратной матрицы
        void swap_colums(int size, double[,] matrix, int c_1, int c_2)
        {
            if ((size <= 0) || (c_1 < 0) || (c_1 >= size) || (c_2 < 0) || (c_2 >= size) || (c_1 == c_2)) return;

            for (int i = 0; i < size; i++) 
            {
                double t = matrix[i, c_1];
                matrix[i, c_1] = matrix[i, c_2];
                matrix[i, c_2] = t;
            }
        }

        // Метод Гаусса с выбором главного элемента по строке 
        // Возврат результата через параметр X_matrix
        int gauss( int size, double[,] A_matrix, double[] B_matrix, double[] X_matrix )
        {
	        if ( size <= 0 ) return 0;

	        // копируем матрицы
	        double[,] a_matrix = new double[size,size];
	        double[]  b_matrix = new double[size];
	        Array.Copy(A_matrix, a_matrix,  size*size);
            Array.Copy(B_matrix, b_matrix,  size);


	        // порядок следования X
	        int[] order = new int[size];
            for (int i = 0; i < size; i++)
                order[i] = i;


	        // формируем треугольную матрицу с перестановкой столбцов
	        for ( int k = 0; k < size; k++ )
	        {

		        int max = k;
		        for ( int i = k + 1; i < size; i++ ) 
		        {
			        if ( Math.Abs( a_matrix[k,i] ) > Math.Abs( a_matrix[k,max] ) ) max = i;
		        }

		        if ( Math.Abs(a_matrix[k,max]) < 0.000001 ) 
		        {
			        return 0;
		        }

		        if (max != k) 
		        {
			        swap_colums( size, a_matrix, k, max );
                    int t = order[k];
                    order[k] = order[max];
                    order[max] = t;
		        }


		        for ( int m = k + 1; m < size; m++ )
		        {
			        double koef = a_matrix[m, k] / a_matrix[k, k];

			        b_matrix[m] = b_matrix[m] - koef * b_matrix[k];

			        a_matrix[m, k] = 0.0;

			        for ( int l = k + 1; l < size; l++ )
			        {
				        a_matrix[m, l] = a_matrix[m, l] - koef * a_matrix[k, l];
			        }			
		        }
				
	        }

	        // расчитываем X с учетом их порядка следования
	        for ( int i = size - 1; i >= 0; i-- )
	        {
		        double sum = 0;
		        for ( int j = i + 1; j < size; j++ )
		        {
			        sum += a_matrix[i, j] * X_matrix[ order[j] ];
		        }

		        X_matrix[ order[i] ] = ( b_matrix[i] - sum ) / a_matrix[i, i];
	        }

	        return 1;
        }

        double[] naim_kvadrat(double[] X, double[] Y, int k) 
        {
            double t;
            int i, j;

            double[,] C = new double[k + 1, k + 1];
            double[] D = new double[k + 1];

            // заполняем матрицу C
            for (int m = 0; m <= 2 * k; m++)
            {
                t = 0;
                for (int p = 0; p < X.Length; p++)
                    t += Math.Pow(X[p], m);

                i = 0;
                j = m;

                while((i<=k)&&(j>=0))
                {
                    if (j <= k) C[i, j] = t;
                    i++;
                    j--;
                }
            }

            // заполняем матрицу D
            for (j = 0; j <= k; j++)
            {
                D[j] = 0;
                for (i = 0; i < X.Length; i++)
                    D[j] += Y[i] * Math.Pow(X[i], j);
            }

            // решаем систему методом Гаусса
            double[] A = new double[k + 1];
            int jop = gauss(k + 1, C, D, A);
            return A;
        }


        private void Form1_Load(object sender, EventArgs e)
        {
            const int N = 5;
            double[] X = new double[N];
            double[] Y = new double[N];

            X[0] = -2;
            X[1] = -1;
            X[2] = 0;
            X[3] = 1;
            X[4] = 2;

            Y[0] = 3;
            Y[1] = 0;
            Y[2] = 10;
            Y[3] = 0;
            Y[4] = 3;


            for (int i = 0; i < N; i++)
                listBox1.Items.Add(X[i]);
            for (int i = 0; i < N; i++)
                listBox2.Items.Add(Y[i]);

            for (int i = 0; i < N; i++)
            { 
                chart1.Series[1].Points.AddXY(X[i], Y[i]);
            }

            double[] A1 = naim_kvadrat(X, Y, 1);
            double[] A2 = naim_kvadrat(X, Y, 2);
            double[] A3 = naim_kvadrat(X, Y, 3);

            for (double i = -2.1; i < 2.1; i += 0.1)
            {
                //double y = newton(X, Y, 0, i);
                double y = lagrange(X, Y, i);
                chart1.Series[0].Points.AddXY(i, y);
                y = calc_polynom(A1, i);
                chart1.Series[2].Points.AddXY(i, y);
                y = calc_polynom(A2, i);
                chart1.Series[3].Points.AddXY(i, y);
                y = calc_polynom(A3, i);
                chart1.Series[4].Points.AddXY(i, y);
            }
            
        }

        private void label1_Click(object sender, EventArgs e)
        {

        }
    }
}
