#ifndef _BACTERIA_H_
#define _BACTERIA_H_
#include "Stokeslet.h"

class Bacteria : public Scene
{
public:
	float time;

	virtual void SceneInit();
};

void Bacteria::SceneInit()
{
	ds = 0.1;
	e = ds / 4;

	showgrid = 0.2;
	xmax = 5;
	ymax = 5;
	zmax = 5;

	//Parameters for the tail
	float tailLength = 2; //Straight length of the tail
	int tailNum = (int)round(tailLength / (ds / 5));
	double omega = -2;   //Rotational angular velocity
	double A = 0.5;       //Max amplitude of oscillator
	double p = 3;         //Periodic

	//Parameters for the sphere
	double radius = 1.0;
	int n1 = (int)round(2 * PI*radius / ds) + 1;
	int n2 = (int)round(PI*radius / ds) + 1;

	num = n1*n2 + tailNum;
	pStokeslet slist = (pStokeslet)malloc(sizeof(Stokeslet)*num);

	//Sphere
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
		{
			double theta = 2 * PI*(i - 1) / n1;
			double phi = PI*(j - 1) / n2;
			slist[i*n2 + j].Position.x = radius*sin(phi)*cos(theta);
			slist[i*n2 + j].Position.y = radius*sin(phi)*sin(theta);
			slist[i*n2 + j].Position.z = radius*cos(phi);
			slist[i*n2 + j].Velocity = Vector(0, 0, 0);
		}
	
	//Tail
	for (int i = 0; i < tailNum; i++)
	{
		double B = A*(exp((double)i / 10) - 1);  //Amplitude
		if (B > A) B = A;
		double x = 0 - radius - ds*i;
		double y = B*sin(p*x);

		//Transform
		double theta = omega*time;
		double vz = y*omega*cos(theta);
		double vy = 0 - y*omega*sin(theta);
		double z = y*sin(theta);
		y = y*cos(theta);
		
		slist[n1*n2 + i].Position = Vector(x, y, z);
		slist[n1*n2 + i].Velocity = Vector(0, vy, vz);
	}

	if (Stokeslets != NULL) free(Stokeslets);
	Stokeslets = slist;
}

#endif