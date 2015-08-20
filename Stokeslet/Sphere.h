#ifndef _SPHERE_H_
#define _SPHERE_H_
#include "Stokeslet.h"

class Sphere :public Scene
{
public:
	virtual void SceneInit();
};

void Sphere::SceneInit()
{
	ds = 0.1;
	e = ds / 4;

	showgrid = 0.1;
	xmax = 4;
	ymax = 4;
	zmax = 4;

	double radius = 1.2;
	int n1 = (int)round(2 * PI*radius / ds) + 1;
	int n2 = (int)round(PI*radius / ds) + 1;
	num = n1*n2;
	pStokeslet slist = (pStokeslet)malloc(sizeof(Stokeslet)*num);
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
		{
			double theta = 2 * PI*(i - 1) / n1;
			double phi = PI*(j - 1) / n2;
			slist[i*n2 + j].Position.x = radius*sin(phi)*cos(theta);
			slist[i*n2 + j].Position.y = radius*sin(phi)*sin(theta);
			slist[i*n2 + j].Position.z = radius*cos(phi);
			slist[i*n2 + j].Velocity = Vector(0.4, 0, 0);
		}
	if (Stokeslets != NULL) free(Stokeslets);
	Stokeslets = slist;
}
#endif