//Level set simulation- 1. Original Function, 2. Reinitialized function without sub-cell fix, 3. Reinitialized function with sub-cell fix
//Domain is the function of [-2,2] x [-2,2]. Grid is 64x64.
// Function 1: phi_0(x,y) = 0.25*x*x+y*y-0.25
// Function 2: phi_0(x,y) = 16*x^4-4*x^2+16*y^4-4*y^2
// Function 3: phi_0 = (x*x+y*y)^2-1

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>

int main(int argc, char* argv[])
{
	//Defining the domain
	double xmin = -2, xmax = 2;
	double ymin = -2, ymax = 2;

	//Defining the number of grid points
	int imax = 64, jmax = 64;

	//Finding out the value of grid spacing
	double dx, dy;
	dx = (xmax - xmin) / (imax - 1);
	dy = (ymax - ymin) / (jmax - 1);

	//Defining the time step
	double dt = 0.5*dx;

	//Calculating the value of X and Y at all the grid points
	std::vector<std::vector<double>> X(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> Y(imax, std::vector<double>(jmax, 0.0));

	for (int j = 0; j < jmax; j++)
	{
		X[0][j] = -2;
		X[imax - 1][j] = 2;
	}
	
	for (int i = 0; i < imax; i++)
	{
		Y[i][0] = -2;
		Y[i][jmax - 1] = 2;
	}

	for (int i = 1; i < imax - 1; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
		X[i][j] = X[i - 1][j] + dx;
		}
	}

	for (int j = 1; j < jmax - 1; j++)
	{
		for (int i = 0; i < imax; i++)
		{
			Y[i][j] = Y[i][j - 1] + dy;
		}
	}

	//Total number of iterations to perform
	int n = 512;

	std::vector<std::vector<double>> phi(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> phi_zero(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> phi_new(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> a(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> b(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> c(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> d(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> a_plus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> b_plus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> c_plus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> d_plus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> a_minus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> b_minus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> c_minus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> d_minus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> s_plus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> s_minus(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> s_zero(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> SDF(imax, std::vector<double>(jmax, 0.0));
	std::vector<std::vector<double>> phi_subfix(imax, std::vector<double>(jmax, 0.0));

	//Initialization of the functions
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			//phi_zero[i][j] = 0.25*(pow(X[i][j], 2) + pow(Y[i][j], 2)) - 0.25;
			phi_zero[i][j] = 16 * (pow(X[i][j], 4)) + 16 * pow(Y[i][j], 4) - 4 * pow(X[i][j], 2) - 4 * pow(Y[i][j], 2);
			//phi_zero[i][j] = pow(X[i][j], 2) + pow(Y[i][j], 2) - 1;
		}
	}

	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			phi[i][j] = phi_zero[i][j];
		}
	}

	//Marching in time
	for (int t = 0; t < n; t++)
	{
		//Reinitialization of phi
		//Calculating the values of a, b, c, and d
		for (int i = 1; i < imax - 1; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				a[i][j] = (phi[i][j] - phi[i - 1][j]) / dx;
				b[i][j] = (phi[i + 1][j] - phi[i][j]) / dx;
				c[i][j] = (phi[i][j] - phi[i][j - 1]) / dy;
				d[i][j] = (phi[i][j + 1] - phi[i][j]) / dy;
			}
		}

		//Calculation of a, b, c and d at boundaries by using Ghost Nodes
		for (int j = 0; j < jmax; j++)
		{
			a[0][j] = (phi[0][j] - phi[imax - 2][j]) / dx;
			b[0][j] = (phi[1][j] - phi[0][j]) / dx;
			a[imax - 1][j] = (phi[imax - 1][j] - phi[imax - 2][j]) / dx;
			b[imax - 1][j] = (phi[1][j] - phi[imax - 1][j]) / dx;
		}
		for (int i = 0; i < imax; i++)
		{
			c[i][0] = (phi[i][0] - phi[i][jmax - 2]) / dy;
			d[i][0] = (phi[i][1] - phi[i][0]) / dy;
			c[i][jmax - 1] = (phi[i][jmax - 1] - phi[i][jmax - 2]) / dy;
			d[i][jmax - 1] = (phi[i][1] - phi[i][jmax - 1]) / dy;
		}
		for (int i = 1; i < imax - 1; i++)
		{
			a[i][0] = (phi[i][0] - phi[i - 1][0]) / dx;
			b[i][0] = (phi[i + 1][0] - phi[i][0]) / dx;
			a[i][jmax - 1] = (phi[i][jmax - 1] - phi[i - 1][jmax - 1]) / dx;
			b[i][jmax - 1] = (phi[i + 1][jmax - 1] - phi[i][jmax - 1]) / dx;
		}
		for (int j = 1; j < jmax - 1; j++)
		{
			c[0][j] = (phi[0][j] - phi[0][j - 1]) / dy;
			c[imax - 1][j] = (phi[imax - 1][j] - phi[imax - 1][j - 1]) / dy;
			d[0][j] = (phi[0][j + 1] - phi[0][j]) / dy;
			d[imax - 1][j] = (phi[imax - 1][j + 1] - phi[imax - 1][j]) / dy;
		}

		//Calculating values of a+, a-, b+, b-, c+, c-, d+, d-, s+, s-
		for (int i = 0; i < imax; i++)
		{
			for (int j = 0; j < jmax; j++)
			{
				a_plus[i][j] = std::max(a[i][j], 0.0);
				a_minus[i][j] = std::min(0.0, a[i][j]);
				b_plus[i][j] = std::max(b[i][j], 0.0);
				b_minus[i][j] = std::min(0.0, b[i][j]);
				c_plus[i][j] = std::max(c[i][j], 0.0);
				c_minus[i][j] = std::min(0.0, c[i][j]);
				d_plus[i][j] = std::max(d[i][j], 0.0);
				d_minus[i][j] = std::min(0.0, d[i][j]);

				if (phi[i][j] > 0.0001)
				{
					s_plus[i][j] = 1;
					s_minus[i][j] = 0;
				}
				else if (phi[i][j] < -0.0001)
				{
					s_plus[i][j] = 0;
					s_minus[i][j] = -1;
				}
				else
				{
					s_plus[i][j] = (phi[i][j]) / (sqrt((pow(phi[i][j], 2)) + pow(dx, 2)));
					s_minus[i][j]= (phi[i][j]) / (sqrt((pow(phi[i][j], 2)) + pow(dx, 2)));
				}

				if (phi[i][j] > 0.0001)
				{
					s_zero[i][j] = 1;
				}
				if (phi[i][j] < -0.0001)
				{
					s_zero[i][j] = -1;
				}
				if (phi[i][j] == 0)
				{
					s_zero[i][j]= (phi[i][j]) / (sqrt((pow(phi[i][j], 2)) + pow(dx, 2)));
				}
			}
		}

		//Finding out the new values of phi
		for (int i = 0; i < imax; i++)
		{
			for (int j = 0; j < jmax; j++)
			{
				phi_new[i][j] = phi[i][j] + dt*(s_plus[i][j] * ((sqrt(std::max(pow(a_plus[i][j], 2), pow(b_minus[i][j], 2))) + std::max(pow(c_plus[i][j], 2), pow(d_minus[i][j], 2))) - 1)) - dt*(s_minus[i][j] * ((sqrt(std::max(pow(a_minus[i][j], 2), pow(b_plus[i][j], 2)) + std::max(pow(c_minus[i][j], 2), pow(d_plus[i][j], 2)))) - 1));
			}
		}


		//Implementing the Sub-Cell Fix
		//Calculating the Signed Distance Function SDF

		for (int i = 1; i < imax - 1; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				if (phi[i + 1][j] * phi[i][j] < 0 || phi[i][j + 1] * phi[i][j] < 0 || phi[i][j] * phi[i - 1][j] < 0 || phi[i][j] * phi[i][j - 1] < 0)
				{
					SDF[i][j] = (2 * dx*phi_zero[i][j]) / (sqrt((pow((phi_zero[i + 1][j] - phi_zero[i - 1][j]), 2)) + (pow((phi_zero[i][j + 1] - phi_zero[i][j - 1]), 2))));
				}
				else
				{
					SDF[i][j] = 0;
				}
			}
		}

		for (int i = 1; i < imax - 1; i++)
		{
			if (phi[i + 1][0] * phi[i][0] < 0 || phi[i][1] * phi[i][0] < 0 || phi[i][0] * phi[i - 1][0] < 0 || phi[i][0] * phi[i][jmax - 2] < 0)
			{
				SDF[i][0] = (2 * dx*phi_zero[i][0]) / (sqrt((pow((phi_zero[i + 1][0] - phi_zero[i - 1][0]), 2)) + pow((phi_zero[i][1] - phi_zero[i][jmax - 2]), 2)));
				phi_subfix[i][0] = phi[i][0] - dt*((s_zero[i][0] * abs(phi[i][0])) - SDF[i][0]) / dx;
			}
			else
			{
				SDF[i][0] = 0;
				phi_subfix[i][0] = phi[i][0] - dt*(s_plus[i][0] * ((sqrt(std::max(pow(a_plus[i][0], 2), pow(b_minus[i][0], 2)) + std::max(pow(c_plus[i][0], 2), pow(d_minus[i][0], 2))) - 1))) - dt*(s_minus[i][0] * ((sqrt(std::max(pow(a_minus[i][0], 2), pow(b_plus[i][0], 2)) + std::max(pow(c_minus[i][0], 2), pow(d_plus[i][0], 2))) - 1)));
			}
			if (phi[i + 1][jmax - 1] * phi[i][jmax - 1] < 0 || phi[i][1] * phi[i][jmax - 1] < 0 || phi[i][jmax - 1] * phi[i - 1][jmax - 1] < 0 || phi[i][jmax - 1] * phi[i][jmax - 2] < 0)
			{
				SDF[i][jmax - 1] = (2 * dx*phi_zero[i][jmax - 1]) / (sqrt((pow((phi_zero[i + 1][jmax - 1] - phi_zero[i - 1][jmax - 1]), 2)) + (pow((phi_zero[i][jmax - 2] - phi_zero[i][1]), 2))));
				phi_subfix[i][jmax - 1] = phi[i][jmax - 1] - dt*((s_zero[i][jmax - 1] * abs(phi[i][jmax - 1])) - SDF[i][jmax - 1]) / dx;
			}
			else
			{
				SDF[i][jmax - 1] = 0;
				phi_subfix[i][jmax - 1] = phi[i][jmax - 1] - dt*(s_plus[i][jmax - 1] * ((sqrt(std::max(pow(a_plus[i][jmax - 1], 2), pow(b_minus[i][jmax - 1], 2)) + std::max(pow(c_plus[i][jmax - 1], 2), pow(d_minus[i][jmax - 1], 2))) - 1))) - dt*(s_minus[i][jmax - 1] * ((sqrt(std::max(pow(a_minus[i][jmax - 1], 2), pow(b_plus[i][jmax - 1], 2)) + std::max(pow(c_minus[i][jmax - 1], 2), pow(d_plus[i][jmax - 1], 2))) - 1)));
			}
		}

		for (int j = 1; j < jmax - 1; j++)
		{
			if (phi[1][j] * phi[0][j] < 0 || phi[0][j + 1] * phi[0][j] < 0 || phi[0][j] * phi[imax - 2][j] < 0 || phi[0][j] * phi[0][j - 1] < 0)
			{
				SDF[0][j] = (2 * dx*phi_zero[0][j]) / (sqrt((pow((phi_zero[1][j] - phi_zero[imax - 2][j]), 2)) + (pow((phi_zero[1][j + 1] - phi_zero[0][j - 1]), 2))));
				phi_subfix[0][j] = phi[0][j] - dt*((s_zero[0][j] * abs(phi[0][j])) - SDF[0][j]) / dx;
			}
			else
			{
				SDF[0][j] = 0;
				phi_subfix[0][j] = phi[0][j] - dt*(s_plus[0][j] * ((sqrt(std::max(pow(a_plus[0][j], 2), pow(b_minus[0][j], 2)) + std::max(pow(c_plus[0][j], 2), pow(d_minus[0][j], 2)))) - 1)) - dt*(s_minus[0][j] * ((sqrt(std::max(pow(a_minus[0][j], 2), pow(b_plus[0][j], 2)) + std::max(pow(c_minus[0][j], 2), pow(d_plus[0][j], 2)))) - 1));
			}
			if (phi[1][j] * phi[imax - 1][j] < 0 || phi[imax - 1][j + 1] * phi[imax - 1][j] < 0 || phi[imax - 1][j] * phi[imax - 2][j] < 0 || phi[imax - 1][j] * phi[imax - 1][j - 1] < 0)
			{
				SDF[imax - 1][j] = (2 * dx*phi_zero[imax - 1][j]) / (sqrt((pow((phi_zero[1][j] - phi_zero[imax - 2][j]), 2) + pow((phi_zero[imax - 1][j + 1] - phi_zero[imax - 1][j - 1]), 2))));
				phi_subfix[imax - 1][j] = phi[imax - 1][j] - dt*((s_zero[imax - 1][j] * abs(phi[imax - 1][j])) - SDF[imax - 1][j]) / dx;
			}
			else
			{
				SDF[imax - 1][j] = 0;
				phi_subfix[imax - 1][j] = phi[imax - 1][j] - dt*(s_plus[imax - 1][j] * ((sqrt(std::max(pow(a_plus[imax - 1][j], 2), pow(b_minus[imax - 1][j], 2)) + std::max(pow(c_plus[imax - 1][j], 2), pow(d_minus[imax - 1][j], 2))))-1)) - dt*(s_minus[imax - 1][j] * ((sqrt(std::max(pow(a_minus[imax - 1][j], 2), pow(b_plus[imax - 1][j], 2)) + std::max(pow(c_minus[imax - 1][j], 2), pow(d_plus[imax - 1][j], 2))))-1));
			}
		}

		//Calulating the value of phi after sub-cell fix by marching in time
		for (int i = 1; i < imax - 1; i++)
		{
			for (int j = 0; j < jmax - 1; j++)
			{
				if (phi[i + 1][j] * phi[i][j] < 0 || phi[i][j + 1] * phi[i][j] < 0 || phi[i][j] * phi[i - 1][j] < 0 || phi[i][j] * phi[i][j - 1] < 0)
				{
					phi_subfix[i][j] = phi[i][j] - dt*((s_zero[i][j] * abs(phi[i][j])) - SDF[i][j]) / dx;
				}
				else
				{
					phi_subfix[i][j] = phi[i][j] - dt*(s_plus[i][j] * ((sqrt(std::max(pow(a_plus[i][j], 2), pow(b_minus[i][j], 2)) + std::max(pow(c_plus[i][j], 2), pow(d_minus[i][j], 2))) - 1))) - dt*(s_minus[i][j] * ((sqrt(std::max(pow(a_minus[i][j], 2), pow(b_plus[i][j], 2)) + std::max(pow(c_minus[i][j], 2), pow(d_plus[i][j], 2))) - 1)));
				}
			}
		}

		phi_subfix[0][0] = phi_subfix[1][0];
		phi_subfix[imax - 1][jmax - 1] = phi_subfix[imax - 2][jmax - 1];
		phi_subfix[0][jmax - 1] = phi_subfix[1][jmax - 1];
		phi_subfix[imax - 1][0] = phi_subfix[imax - 2][0];

		//Updating the phi
		//Un-comment this part of code when sub-cell fix required
		/*for (int i = 0; i < imax; i++)
		{
			for (int j = 0; j < jmax; j++)
			{
				phi[i][j] = phi_subfix[i][j];
			}
		}*/


	}

	//Writing the files for plots
	std::ofstream write_X("X.dat");
	write_X.setf(std::ios::scientific);
	write_X.precision(16);
	assert(write_X.is_open());
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			write_X << X[i][j] << " ";
		}
		write_X << "\n";
	}
	write_X.close();

	std::ofstream write_Y("Y.dat");
	write_Y.setf(std::ios::scientific);
	write_Y.precision(16);
	assert(write_Y.is_open());
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			write_Y << Y[i][j] << " ";
		}
		write_Y << "\n";
	}
	write_Y.close();

	std::ofstream write_phi_zero("phi_zero.dat");
	write_phi_zero.setf(std::ios::scientific);
	write_phi_zero.precision(16);
	assert(write_phi_zero.is_open());
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			write_phi_zero << phi_zero[i][j] << " ";
		}
		write_phi_zero << "\n";
	}
	write_phi_zero.close();

	std::ofstream write_phi("phi.dat");
	write_phi.setf(std::ios::scientific);
	write_phi.precision(16);
	assert(write_phi.is_open());
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			write_phi << phi[i][j] << " ";
		}
		write_phi << "\n";
	}
	write_phi.close();

	return 0;
}