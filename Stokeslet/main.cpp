#include <iostream>
#include <string>
#include "Stokeslet.h"
#include "Bacteria.h"
using namespace std;

int main()
{
	Bacteria *scene = new Bacteria();
	scene->time = 2;
	scene->SceneInit();
	scene->Solve();
	scene->Output("stokeslets.out", "flow.out");
	delete scene;

	return 0;
}