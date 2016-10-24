#ifndef FLUID_SOLVER_2D_H_
#define FLUID_SOLVER_2D_H_

#include <string>

class FluidSolver2D {
	typedef int** Mat2Di;
	typedef float** Mat2Df;

private:
	//----------------------------------------------------------------------
	// Grid Attributes
	//----------------------------------------------------------------------
	// nx
	int gridWidth;
	// ny
	int gridHeight;
	// distance between each grid cell
	float dx;
	// grid of cell labels, size (nx, ny)
	Mat2Di label;
	// grid of pressures, size (nx, ny)
	Mat2Df p;
	// grid of vel x component, size (nx+1, ny)
	Mat2Df u;
	// grid of vel y component, size (nx, ny+1)
	Mat2Df v;

	//----------------------------------------------------------------------
	// Simulation Attributes
	//----------------------------------------------------------------------
	float dt;

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	/*
	Builds initial grid of dimensions (width, height) that contains the initial
	geometry for the system to simulate. It reads the initial geometry from
	the specified input file parameter.
	Args:
	width/height - grid dimensions
	geomFile - the file containing the geometry
	grid - the 2D array to put the initial grid in
	*/
	void readInGeom(int, int, std::string, Mat2Di);


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

	/*
	Initializes the solver by reading in and constructing initial
	grid based on the given initial geometry file.
	Args:
	initialGemoetryFile - name of the .txt file containing initial geometry
	*/
	void init(std::string);

};

#endif //FLUID_SOLVER_2D_H_