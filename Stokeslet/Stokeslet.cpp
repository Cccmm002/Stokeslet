#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "Stokeslet.h"
#include "mgmres.h"
using namespace std;

void Scene::Solve()
{
	double e2 = e*e;

	//Calculate the matrix M.
	double *M = (double*)malloc(sizeof(double)*num*num * 9);
	int *iM = (int*)malloc(sizeof(int)*num*num * 9);
	int *jM = (int*)malloc(sizeof(int)*num*num * 9);
	for (int b = 0; b < num; b++)
		for (int k = 0; k < num; k++)
		{
			Vector r = Stokeslets[b].Position - Stokeslets[k].Position;
			double r2 = r.AbsSquare();
			double H1 = (double)1 / (8 * PI*sqrt(r2 + e2)) + e2 / (8 * PI*pow(r2 + e2, (double)3 / 2));
			double H2 = (double)1 / (8 * PI*pow(r2 + e2, (double)3 / 2));

			double M1 = H1;
			double M2xx = H2*r.x*r.x;
			double M2yy = H2*r.y*r.y;
			double M2zz = H2*r.z*r.z;
			double M2xy = H2*r.x*r.y;
			double M2yz = H2*r.y*r.z;
			double M2xz = H2*r.x*r.z;

			//The third-party GMRES module is designed for sparse matrix and the input matrix is in
			//sparse triplet form, which is using array M to store the matrix, array iM to store the
			//row of each entry in the matrix and jM to store the column of each entry. For details,
			//please read the description in the function "mgmres_st" in "mgmres.cpp".
			//Using O3 optimization of your compiler can make the following block of code run a
			//little bit faster.

			//M=[M1+M2xx M2xy M2xz; M2xy M1+M2yy M2yz; M2xz M2yz M1+M2zz]
			M[b*(3 * num) + k] = M1 + M2xx;
			iM[b*(3 * num) + k] = b; jM[b*(3 * num) + k] = k;
			M[b*(3 * num) + k + num] = M2xy;
			iM[b*(3 * num) + k + num] = b; jM[b*(3 * num) + k + num] = k + num;
			M[b*(3 * num) + k + num + num] = M2xz;
			iM[b*(3 * num) + k + num + num] = b; jM[b*(3 * num) + k + num + num] = k + num + num;
			M[(b + num)*(3 * num) + k] = M2xy;
			iM[(b + num)*(3 * num) + k] = b + num; jM[(b + num)*(3 * num) + k] = k;
			M[(b + num)*(3 * num) + k + num] = M1 + M2yy;
			iM[(b + num)*(3 * num) + k + num] = b + num; jM[(b + num)*(3 * num) + k + num] = k + num;
			M[(b + num)*(3 * num) + k + num + num] = M2yz;
			iM[(b + num)*(3 * num) + k + num + num] = b + num; jM[(b + num)*(3 * num) + k + num + num] = k + num + num;
			M[(b + num + num)*(3 * num) + k] = M2xz;
			iM[(b + num + num)*(3 * num) + k] = b + num + num; jM[(b + num + num)*(3 * num) + k] = k;
			M[(b + num + num)*(3 * num) + k + num] = M2yz;
			iM[(b + num + num)*(3 * num) + k + num] = b + num + num; jM[(b + num + num)*(3 * num) + k + num] = k + num;
			M[(b + num + num)*(3 * num) + k + num + num] = M1 + M2zz;
			iM[(b + num + num)*(3 * num) + k + num + num] = b + num + num; jM[(b + num + num)*(3 * num) + k + num + num] = k + num + num;
		}
	//Initialize array U and F.
	double *U = (double*)malloc(sizeof(double)*num * 3);
	double *F = (double*)malloc(sizeof(double)*num * 3);
	memset(F, 0, sizeof(double)*num * 3);   //Initial 0 guess
	for (int i = 0; i < num; i++)
	{
		U[i] = Stokeslets[i].Velocity.x;
		U[i + num] = Stokeslets[i].Velocity.y;
		U[i + num + num] = Stokeslets[i].Velocity.z;
	}

	//Using GMRES to solve the equation.
	mgmres_st(num * 3, 9 * num*num, iM, jM, M, F, U, 30, num * 3 - 1, 1e-8, 1e-8);

	//Obtain the solving result
	for (int i = 0; i < num; i++)
		Stokeslets[i].Force = Vector(F[i], F[i + num], F[i + num + num]);

	//Free memory
	free(F); free(U);
	free(M); free(iM); free(jM);


	//Prepare to calculate flow velocity of points in flow field.
	xgrid = (int)floor(2 * (float)xmax / showgrid);
	ygrid = (int)floor(2 * (float)ymax / showgrid);
	zgrid = (int)floor(2 * (float)zmax / showgrid);
	if (points != NULL) free(points);
	points = (pFlowPoint)malloc(sizeof(FlowPoint)*xgrid*ygrid*zgrid);

	//In order to use OpenMP to accelerate, follow the instruction of your compiler to turn on OpenMP.
#pragma omp parallel for schedule(dynamic)
	for (int a = 0; a < xgrid; a++)
		for (int b = 0; b < ygrid; b++)
			for (int c = 0; c < zgrid; c++)
			{
				pFlowPoint curPoint = points + a*ygrid*zgrid + b*zgrid + c;
				Vector position = Vector((0 - xmax) + showgrid*a, (0 - ymax) + showgrid*b, (0 - zmax) + showgrid*c);
				curPoint->Position = position;
				curPoint->Velocity = Vector(0, 0, 0);

				for (int k = 0; k < num; k++)
				{
					Vector r = position - Stokeslets[k].Position;
					double r2 = r.AbsSquare();
					double H1 = (double)1 / (8 * PI*sqrt(r2 + e2)) + e2 / (8 * PI*pow(r2 + e2, (double)3 / 2));
					double H2 = (double)1 / (8 * PI*pow(r2 + e2, (double)3 / 2));
					double fk = Stokeslets[k].Force*r;
					curPoint->Velocity = curPoint->Velocity + H1*Stokeslets[k].Force + H2*(fk*r);
				}
			}
}

void Scene::Output(string stokeslets_path, string flow_path)
{
	if (Stokeslets != NULL && stokeslets_path != "")
	{
		//Output: Stokeslets
		//=====================================
		//First line: num
		//Then: x, y, z coordinates of Stokeslets (lines 2-4);
		//The following 3 lines: Velocity of Stokeslets(lines 5-7);
		//Then forces of Stokeslets(lines 8-10);
		ofstream file_s;
		file_s.open(stokeslets_path);
		file_s << num << endl;

		ostringstream posx, posy, posz;
		ostringstream speedx, speedy, speedz;
		ostringstream forcex, forcey, forcez;
		for (int i = 0; i < num; i++)
		{
			if (i != 0)
			{
				posx << ','; posy << ','; posz << ',';
				speedx << ','; speedy << ','; speedz << ',';
				forcex << ','; forcey << ','; forcez << ',';
			}
			posx << Stokeslets[i].Position.x; posy << Stokeslets[i].Position.y; posz << Stokeslets[i].Position.z;
			speedx << Stokeslets[i].Velocity.x; speedy << Stokeslets[i].Velocity.y; speedz << Stokeslets[i].Velocity.z;
			forcex << Stokeslets[i].Force.x; forcey << Stokeslets[i].Force.y; forcez << Stokeslets[i].Force.z;
		}
		file_s << posx.str() << endl << posy.str() << endl << posz.str() << endl;
		file_s << speedx.str() << endl << speedy.str() << endl << speedz.str() << endl;
		file_s << forcex.str() << endl << forcey.str() << endl << forcez.str() << endl;

		file_s.close();
	}

	if (points != NULL && flow_path != "")
	{
		//Output: Flow velocity
		//=====================================
		//First line: xgrid, zgrid, zgrid.
		//All arranged by z, first position of points, then velocity.
		//For example, if it is xgrid*ygrid*zgrid, then for every z, there
		//is an xgrid*ygrid matrix. There are zgrid matrics.
		//Totally, there are (2*3*zgrid + 1) lines in the file, with no blank lines.
		ofstream file_f;
		file_f.open(flow_path);
		file_f << xgrid << ',' << ygrid << ',' << zgrid << endl;

		ostringstream *fpx = new ostringstream[zgrid], *fpy = new ostringstream[zgrid], *fpz = new ostringstream[zgrid];
		ostringstream *fvx = new ostringstream[zgrid], *fvy = new ostringstream[zgrid], *fvz = new ostringstream[zgrid];
		for (int y = 0; y < ygrid; y++)
		{
			for (int z = 0; z < zgrid; z++)
				for (int x = 0; x < xgrid; x++)
				{
					pFlowPoint curPoint = points + x*ygrid*zgrid + y*zgrid + z;
					if (x != 0)
					{
						fpx[z] << ','; fpy[z] << ','; fpz[z] << ',';
						fvx[z] << ','; fvy[z] << ','; fvz[z] << ',';
					}
					fpx[z] << curPoint->Position.x; fpy[z] << curPoint->Position.y; fpz[z] << curPoint->Position.z;
					fvx[z] << curPoint->Velocity.x; fvy[z] << curPoint->Velocity.y; fvz[z] << curPoint->Velocity.z;
				}
			for (int z = 0; z < zgrid; z++)
			{
				fpx[z] << endl; fpy[z] << endl; fpz[z] << endl;
				fvx[z] << endl; fvy[z] << endl; fvz[z] << endl;
			}
		}

		for (int i = 0; i < zgrid; i++) file_f << fpx[i].str();
		for (int i = 0; i < zgrid; i++) file_f << fpy[i].str();
		for (int i = 0; i < zgrid; i++) file_f << fpz[i].str();
		for (int i = 0; i < zgrid; i++) file_f << fvx[i].str();
		for (int i = 0; i < zgrid; i++) file_f << fvy[i].str();
		for (int i = 0; i < zgrid; i++) file_f << fvz[i].str();

		/*delete fpx; delete fpy; delete fpz;
		delete fvx; delete fvy; delete fvz;*/

		file_f.close();
	}
}
