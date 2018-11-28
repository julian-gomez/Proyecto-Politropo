//En este codigo pretendemos generar los datos del indice M_anisotropica/M_isotropica usando el codigo anisotropico y teniendo en cuenta que el caso sin anisotropia se reduce haciendo alpha = beta y m = n. Similar al anterior, pero en este fijaremos n y variaremos m entre 0.1 y 4.5
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//Se resolvera usando el metodo de Runge Kutta 4 orden
double h = 0.00001; //Paso de iteracion
double func1(double z, double w1, double w2, double m, double n, double cr, double ct); //cr constante radial; ct constante tangencial
double func2(double z, double w1, double w2, double m, double n, double cr, double ct);
//Condiciones Iniciales
double z_0 = 0.00001;
double w_0 = 1.0;
double der_0 = 0.0;

//Constantes

//Politropo Radial
double beta = 1.0;
double longitud = 30.0;
double delta = 0.01;
int numPuntos = int(longitud/delta);
//Politropo Tangencial

//Definimos hasta donde correremos el indice n
double n = 3.0;

double alpha = 1.0;

//Densidad central
double p_c = 1.0;
int main ()
{
	//Sistema 2 hace referencia al sistema donde alpha = beta y n = m
	int i,j,k;
	double m = 0.5; //m inicial
	double indice; //Relacion de masa anisotropica sobre masa isotropica
	int req1 = 0;
	int req2 = 0;
	double z = z_0;
	double z2 = z_0;
	double w = w_0;
	double w_2 = w_0;
	double resta;
	double w_prime = der_0;
	double w_prime2 = der_0;
	double k1_1,k2_1,k3_1,k4_1,k1_2,k2_2,k3_2,k4_2,promedio1,promedio2;
	ofstream archivo_salida("masas.txt");
	
	for(i = 0; i < numPuntos-1; i++)
	{
		z = z_0;
		w = w_0;
		w_prime = der_0;
		z2 = z_0;
		w_2 = w_0;
		w_prime2 = der_0;
		resta = -1.0*beta*(m+1.0)*pow(p_c,1.0/m); //Phi_c - F_c
		while (req2 == 0)
		{
			k1_1 = func1(z2,w_2,w_prime2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);
			k1_2 = func2(z2,w_2,w_prime2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);
	
			k2_1 = func1(z2 + h/2.0, w_2 + h/2.0*k1_1, w_prime2 + h/2.0*k1_2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);
			k2_2 = func2(z2 + h/2.0, w_2 + h/2.0*k1_1, w_prime2 + h/2.0*k1_2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);

			k3_1 = func1(z2 + h/2.0, w_2 + h/2.0*k2_1, w_prime2 + h/2.0*k2_2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);
			k3_2 = func2(z2 + h/2.0, w_2 + h/2.0*k2_1, w_prime2 + h/2.0*k2_2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);

			k4_1 = func1(z2 + h, w_2 + h*k3_1, w_prime2 + h*k3_2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);
			k4_2 = func2(z2 + h, w_2 + h*k3_1, w_prime2 + h*k3_2,(i+1)*delta + m,(i+1)*delta + m, beta, beta);

			promedio1 = (k1_1 + 2.0*k2_1 + 2.0*k3_1 + k4_1)/6.0;
			promedio2 = (k1_2 + 2.0*k2_2 + 2.0*k3_2 + k4_2)/6.0;

			z2 = z2 + h;
			w_2 = w_2 + h*promedio1;
			w_prime2 = w_prime2 + h*promedio2;
		
			if (w_2 < 0.0001)
			{
				req2 = 1;
			} 
		}
		while (req1 == 0)
		{
			k1_1 = func1(z,w,w_prime,(i+1)*delta + m,n, beta, alpha);
			k1_2 = func2(z,w,w_prime,(i+1)*delta + m,n, beta, alpha);
	
			k2_1 = func1(z + h/2.0, w + h/2.0*k1_1, w_prime + h/2.0*k1_2,(i+1)*delta + m,n, beta, alpha);
			k2_2 = func2(z + h/2.0, w + h/2.0*k1_1, w_prime + h/2.0*k1_2,(i+1)*delta + m,n, beta, alpha);

			k3_1 = func1(z + h/2.0, w + h/2.0*k2_1, w_prime + h/2.0*k2_2,(i+1)*delta + m,n, beta, alpha);
			k3_2 = func2(z + h/2.0, w + h/2.0*k2_1, w_prime + h/2.0*k2_2,(i+1)*delta + m,n, beta, alpha);

			k4_1 = func1(z + h, w + h*k3_1, w_prime + h*k3_2,(i+1)*delta + m,n, beta, alpha);
			k4_2 = func2(z + h, w + h*k3_1, w_prime + h*k3_2,(i+1)*delta + m,n, beta, alpha);

			promedio1 = (k1_1 + 2.0*k2_1 + 2.0*k3_1 + k4_1)/6.0;
			promedio2 = (k1_2 + 2.0*k2_2 + 2.0*k3_2 + k4_2)/6.0;

			z = z + h;
			w = w + h*promedio1;
			w_prime = w_prime + h*promedio2;
		
			if (w < 0.0001)
			{
				req1 = 1;
			} 
		}
		cout << z << " " << z2 << endl;
		
		indice = (pow(z,2.0)*(w_prime + 2.0/(resta*z)*(alpha*pow(p_c,1.0/n)*pow(w,((i+1)*delta + m)/n) - beta*pow(p_c,1.0/((i+1)*delta + m))*w)))/(pow(z2,2.0)*w_prime2);
		archivo_salida << (i+1)*delta + m << " " << indice << endl;
		req1 = 0.0;
		req2 = 0.0;
	}
	return 0;
}

double func1(double z, double w1, double w2, double m, double n, double cr, double ct)
{
	return w2;
}

double func2(double z, double w1, double w2, double m, double n, double cr, double ct)
{
	double diferencia = -1.0*cr*(m+1.0)*pow(p_c,1.0/m); //Phi_c - F_c
	return -1.0*pow(w1,m)-2.0*w2/z + 2.0/diferencia*pow(z,-2.0)*(ct*pow(p_c,1.0/n)*pow(w1,m/n) - cr*w1*pow(p_c,1.0/m)) + 2.0/diferencia*pow(z,-2.0)*(z*m*w2*(ct/n*pow(p_c,1.0/n)*pow(w1,m/n - 1.0) - cr/m*pow(p_c,1.0/m)));
}

