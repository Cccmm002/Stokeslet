#ifndef _STOKESLET_H_
#define _STOKESLET_H_

#include <cmath>
#include <string>
using namespace std;

#define NULL 0
#define PI 3.1415926535897932384626433832795

typedef struct Vector
{
	double x, y, z;

	Vector(double a = 0, double b = 0, double c = 0)
	{
		x = a; y = b; z = c;
	}

	inline double AbsSquare()
	{
		return x*x + y*y + z*z;
	}
}*pVector;

typedef struct Stokeslet
{
	Vector Position;
	Vector Velocity;
	Vector Force;
}*pStokeslet;

typedef struct FlowPoint
{
	Vector Position;
	Vector Velocity;
}*pFlowPoint;

class Scene
{
public:
	//Constants and parameters
	double ds; //Stokeslet seperation
	double e;  //Regularization parameter

	double showgrid; //Interval between two points in the flow field
	double xmax, ymax, zmax; //Flow scene size: from -xmax to xmax; -ymax to ymax; -zmax to zmax
	int xgrid, ygrid, zgrid; //Number of points in each direction

	int num;  //Number of Stokeslets in the scene
	pStokeslet Stokeslets; //List of Stokeslets
	pFlowPoint points; //List of points selected in flow field

	//Methods
	virtual void SceneInit() = 0;
	void Solve();
	void Output(string stokeslets_path, string flow_path);

	//Constructor
	Scene() 
	{ 
		Stokeslets = NULL;
		points = NULL;
	}

	~Scene()
	{
		if (points != NULL)
		{ free(points); points = NULL; }
		if (Stokeslets != NULL)
		{ free(Stokeslets); Stokeslets = NULL; }
	}
};


//Operator overloading

inline Vector operator +(const Vector &a, const Vector &b)
{
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector operator -(const Vector &a, const Vector &b)
{
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vector operator *(const double &a, const Vector &b)
{
	return Vector(a*b.x, a*b.y, a*b.z);
}

inline Vector operator *(const int &a, const Vector &b)
{
	return Vector(a*b.x, a*b.y, a*b.z);
}

inline double operator *(const Vector &a, const Vector &b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
#endif