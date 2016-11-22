#ifndef FLUID_SOLVER_2D_H_
#define FLUID_SOLVER_2D_H_

#include <string>
#include <vector>
#include <fstream>
#include <climits>

#include "SimUtil.h"


class FluidSolver2D {

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
	SimUtil::Mat2Di m_label;
	// pressure and velocity are held in a MAC grid so that
	// p(i, j, k) = p_i_j_k
	// u(i, j, k) = u_i-1/2_j_k
	// v(i, j, k) = v_i_j-1/2_k
	// grid of pressures, size (nx, ny)
	SimUtil::Mat2Df m_p;
	// grid of vel x component, size (nx+1, ny)
	SimUtil::Mat2Df m_u;
	// grid of vel y component, size (nx, ny+1)
	SimUtil::Mat2Df m_v;

	// saved grid of vel x component for FLIP update, size (nx+1, ny)
	SimUtil::Mat2Df m_uSaved;
	// saved grid of vel y component for FLIP update, size (nx, ny+1)
	SimUtil::Mat2Df m_vSaved;

	//----------------------------------------------------------------------
	// Simulation Attributes
	//----------------------------------------------------------------------

	const int VEL_UNKNOWN = INT_MIN;
	// number of particles to seed in each cell at start of sim
	const int PARTICLES_PER_CELL = 4;
	// acceleration due to gravity
	const SimUtil::Vec2 GRAVITY = { 0.0f, -9.81f };
	// simulation time step
	float m_dt;

	//----------------------------------------------------------------------
	// Particle-related Members
	//----------------------------------------------------------------------

	// list of all particles in the simulation
	std::vector<SimUtil::Particle2D> *m_particles;

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	// solver steps

	void seedParticles(int, std::vector<SimUtil::Particle2D>*);
	void labelGrid();
	void particlesToGrid();
	void extrapolateGridFluidData(SimUtil::Mat2Df, int, int, int);
	void saveVelocityGrids();
	void applyBodyForces();

	// helper functions
	void initVelGrid(SimUtil::Mat2Df, int, int, float);
	double trilinearHatKernel(SimUtil::Vec2);
	double hatFunction(double);
	double quadBSplineKernel(SimUtil::Vec2);
	bool isFluid(int, int);
	bool checkNeighbors(SimUtil::Mat2Di, int[2], int[2], int[][2], int, int);

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
	grid based on the given initial geometry file. The solver will save particle
	data at each time step to the given output file.
	Args:
	initialGemoetryFile - name of the .txt file containing initial geometry
	*/
	void init(std::string);

	/*
	Steps the simulation forward dt.
	*/
	void step();

	/*
	Saves all particles currently in the simulation using the given file stream
	Outputs in csv format where each particle position is an entry in the row:
	"0.234 0.154, ... , \n"
	Args:
	particleOut - pointer to file stream to use for output
	*/
	void saveParticleData(std::ofstream*);
};

#endif //FLUID_SOLVER_2D_H_