#ifndef FLUID_SOLVER_2D_H_
#define FLUID_SOLVER_2D_H_

#include <string>

class FluidSolver2D {
	typedef int** Mat2Di;
	typedef float** Mat2Df;

private:
	


public:
	/*
	Creates a new 2D fluid solver.
	Args:
	width - width of the grid to use
	height - height of the grid to use
	dx - the grid cell width
	dt - the timestep to use
	*/
	FluidSolver2D(int, int, float, float);
	~FluidSolver2D();

	void init(std::string initialGeometryFile);

};

#endif //FLUID_SOLVER_2D_H_