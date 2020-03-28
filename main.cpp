#include <iostream>
#include <math.h>

using namespace std;

// ������� -x^2 + 7x -4lnx -7


// ��������� 1 ������ �� [0;1]
double koren1(double e) {
	double xn = 0, xk;
	int i = 0;

	do {
		xk = xn;
		xn = exp((xk*xk-7*xk+7)/(-4));
		i++;
	}while(abs(xn-xk) >= 0.031*e);

	cout << "ITER " << i << endl;

	return xn;
}

// ��������� 2 ������ �� [1;2]
double koren2(double e) {
	double xn = 1.5, xk;
	int i = 0;

	do {
		xk = xn;
		xn = (xk*xk+4*log(xk)+7)/7;
		i++;
	}while(abs(xn-xk) >= 0.16*e);

	cout << "ITER " << i << endl;

	return xn;
}

// ��������� 3 ������ �� [3;4]
double koren3(double e) {
	double xn = 3.5, xk;
	int i = 0;

	do {
		xk = xn;
		xn = exp((xk*xk-7*xk+7)/(-4));
		i++;
	}while(abs(xn-xk) >= 0.25*e);

	cout << "ITER " << i << endl;

	return xn;
}


int main() {
	
	cout << koren2(0.01);
	cin.get();
}