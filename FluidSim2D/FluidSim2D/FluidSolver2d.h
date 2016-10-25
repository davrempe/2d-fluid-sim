#ifndef FLUID_SOLVER_2D_H_
#define FLUID_SOLVER_2D_H_

#include <string>
#include <vector>

#include "SimUtil.h"

class FluidSolver2D {
	typedef int** Mat2Di;
	typedef float** Mat2Df;

private:
	//----------------------------------------------------------------------
	// Grid Attributes
	//----------------------------------------------------------------------

	// nx
	int m_gridWidth;
	// ny
	int m_gridHeight;
	// distance between each grid cell
	float m_dx;
	// grid of cell labels, size (nx, ny)
	Mat2Di m_label;
	// grid of pressures, size (nx, ny)
	Mat2Df m_p;
	// grid of vel x component, size (nx+1, ny)
	Mat2Df m_u;
	// grid of vel y component, size (nx, ny+1)
	Mat2Df m_v;

	//----------------------------------------------------------------------
	// Simulation Attributes
	//----------------------------------------------------------------------

	// number of particles to seed in each cell at start of sim
	const int PARTICLES_PER_CELL = 4;
	// simulation time step
	float m_dt;

	//----------------------------------------------------------------------
	// Particle-related Members
	//----------------------------------------------------------------------
	std::vector<SimUtil::Particle2D> *m_particles;

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	
	void readInGeom(int, int, std::string, Mat2Di);
	void seedParticles(int, std::vector<SimUtil::Particle2D>*);
	SimUtil::Vec2 getCellLocation(int, int);

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