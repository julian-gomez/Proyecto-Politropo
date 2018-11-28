#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//Se resolvera usando el metodo de Runge Kutta 4 orden
double h = 0.00001; //Paso de iteracion
double func1(double z, double w1, double w2);
double func2(double z, double w1, double w2);
//Condiciones Iniciales
double z_0 = 0.00001;
double w_0 = 1.0;
double der_0 = 0.0;

//Constantes

//Politropo Radial
double m = 10.0;
double beta = 0.8;

//Politropo Tangencial
double n = 10.0;
double alpha = 1.0;

//Densidad central
double p_c = 1.0;
double diferencia = -1.0*beta*(m+1.0)*pow(p_c,1.0/m); //Phi_c - F_c
int main ()
{
	int req = 0;
	double z = z_0;
	double w = w_0;
	double w_prime = der_0;
	double k1_1,k2_1,k3_1,k4_1,k1_2,k2_2,k3_2,k4_2,promedio1,promedio2;
	ofstream archivo_salida("casoAnisotropo.txt");
	
	archivo_salida << z << " " << w << " " << w_prime << endl;

	while (req == 0)
	{
		k1_1 = func1(z,w,w_prime);
		k1_2 = func2(z,w,w_prime);
	
		k2_1 = func1(z + h/2.0, w + h/2.0*k1_1, w_prime + h/2.0*k1_2);
		k2_2 = func2(z + h/2.0, w + h/2.0*k1_1, w_prime + h/2.0*k1_2);

		k3_1 = func1(z + h/2.0, w + h/2.0*k2_1, w_prime + h/2.0*k2_2);
		k3_2 = func2(z + h/2.0, w + h/2.0*k2_1, w_prime + h/2.0*k2_2);

		k4_1 = func1(z + h, w + h*k3_1, w_prime + h*k3_2);
		k4_2 = func2(z + h, w + h*k3_1, w_prime + h*k3_2);

		promedio1 = (k1_1 + 2.0*k2_1 + 2.0*k3_1 + k4_1)/6.0;
		promedio2 = (k1_2 + 2.0*k2_2 + 2.0*k3_2 + k4_2)/6.0;

		z = z + h;
		w = w + h*promedio1;
		w_prime = w_prime + h*promedio2;
		if (w >= 0.5)
		{
			archivo_salida << z << " " << w << " " << w_prime << endl;
		}	
		else if (w < 0.5)
		{
			req = 1;
		} 
	}
	return 0;
}

double func1(double z, double w1, double w2)
{
	return w2;
}

double func2(double z, double w1, double w2)
{
	return -1.0*pow(w1,m)-2.0*w2/z - 2.0/diferencia*pow(z,-2.0)*(alpha*pow(p_c,1.0/n)*pow(w1,m/n) - beta*w1*pow(p_c,1.0/m)) - 2.0/diferencia*pow(z,-2.0)*(z*m*w2*(alpha/n*pow(p_c,1.0/n)*pow(w1,m/n - 1.0) - beta/m*pow(p_c,1.0/m)));
}









