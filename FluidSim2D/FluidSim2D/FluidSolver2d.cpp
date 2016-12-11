#include <iostream>
#include <exception>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>

#include "FluidSolver2d.h"
#include "SimUtil.h"

#define DEBUG 0

using namespace SimUtil;

//----------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------

FluidSolver2D::FluidSolver2D(int width, int height, float dx, float dt){
	m_gridWidth = width;
	m_gridHeight = height;
	m_dx = dx;
	m_dt = dt;

	m_particles = new std::vector<Particle2D>();
}

//----------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------

FluidSolver2D::~FluidSolver2D(){
	// clean up grid
	deleteGrid2D<int>(m_gridWidth, m_gridHeight, m_label);
	deleteGrid2D<float>(m_gridWidth, m_gridHeight, m_p);
	deleteGrid2D<float>(m_gridWidth + 1, m_gridHeight, m_u);
	deleteGrid2D<float>(m_gridWidth + 1, m_gridHeight, m_uSaved);
	deleteGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_v);
	deleteGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_vSaved);

	delete m_particles;
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

void FluidSolver2D::init(std::string initialGeometryFile){
	// set up the grid for simulation by initializing arrays
	m_label = initGrid2D<int>(m_gridWidth, m_gridHeight);
	m_p = initGrid2D<float>(m_gridWidth, m_gridHeight);
	m_u = initGrid2D<float>(m_gridWidth + 1, m_gridHeight);
	m_uSaved = initGrid2D<float>(m_gridWidth + 1, m_gridHeight);
	m_v = initGrid2D<float>(m_gridWidth, m_gridHeight + 1);
	m_vSaved = initGrid2D<float>(m_gridWidth, m_gridHeight + 1);

	// init vel grids with unknown label value
	initGridValues<float>(m_u, m_gridWidth + 1, m_gridHeight, VEL_UNKNOWN);
	initGridValues<float>(m_v, m_gridWidth, m_gridHeight + 1, VEL_UNKNOWN);
	if (DEBUG) printGrid2D<float>(m_gridWidth + 1, m_gridHeight, m_u);
	if (DEBUG) printGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_v);

	// read in initial geometry to populate label grid
	readInGeom(m_gridWidth, m_gridHeight, initialGeometryFile, m_label);
	if (DEBUG) printGrid2D<int>(m_gridWidth, m_gridHeight, m_label);

	// seed particles using label grid
	seedParticles(PARTICLES_PER_CELL, m_particles);
}

void FluidSolver2D::step() {
	// update the grid labels
	labelGrid();
	if (DEBUG) printGrid2D<int>(m_gridWidth, m_gridHeight, m_label);
	// transfer particle vel to grid
	particlesToGrid();
	if (DEBUG) printGrid2D<float>(m_gridWidth + 1, m_gridHeight, m_u);
	if (DEBUG) printGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_v);
	// extrapolate fluid data out one cell for accurate divergence calculations
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, 2);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, 2);
	if (DEBUG) printGrid2D<float>(m_gridWidth + 1, m_gridHeight, m_u);
	if (DEBUG) printGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_v);
	// save copy of current grid velocities for FLIP update
	saveVelocityGrids();
	// apply body forces on grid (gravity)
	applyBodyForces();
	if (DEBUG) printGrid2D<float>(m_gridWidth + 1, m_gridHeight, m_u);
	if (DEBUG) printGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_v);
	// solve for pressure
	pressureSolve();
	if (DEBUG) printGrid2D<float>(m_gridWidth, m_gridHeight, m_p);
	// apply pressure force
	applyPressure();
	if (DEBUG) printGrid2D<float>(m_gridWidth + 1, m_gridHeight, m_u);
	if (DEBUG) printGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_v);
	// TODO is grid getting cleared before this?
	// transfer grid velocities back to particles
	gridToParticles(PIC_WEIGHT);
	// advect particles
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, m_gridWidth);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, m_gridHeight);
	advectParticles(ADVECT_MAX);
	// detect particles that have penetrated solid boundary and move back inside fluid
	cleanupParticles(m_dx / 4.0f);
}

void FluidSolver2D::saveParticleData(std::ofstream *particleOut) {
	if (particleOut->is_open()) {
		// print out all particle data on same line, each pos separated by ","
		size_t numParticles = m_particles->size();
		if (numParticles > 0) {
			for (int i = 0; i < numParticles - 1; i++) {
				(*particleOut) << m_particles->at(i).pos.x << " " << m_particles->at(i).pos.y << ",";
			}
			(*particleOut) << m_particles->at(numParticles - 1).pos.x << " " << m_particles->at(numParticles - 1).pos.y << "\n";
		} else {
			(*particleOut) << "\n";
		}
	}
}


//----------------------------------------------------------------------
// Private Main Solver Step Functions
//----------------------------------------------------------------------

/*
Seeds the initial simulation particles. Particles are created for each fluid-labeled
cell in a random-jittered subgrid pattern.
Args:
particlesPerCell - number of particles to seed in each fluid cell.
particleList - list to place the new particles in
*/
void FluidSolver2D::seedParticles(int particlesPerCell, std::vector<Particle2D> *particleList) {
	// set random seed
	srand(time(NULL));
	// go through all cells marked fluid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label[i][j] == SimUtil::FLUID) {
				// seed randomly in 2x2 subgrid of the cell
				Vec2 cellCenter = getGridCellPosition(i, j, m_dx);
				Vec2 subCenters[] = {
					Vec2(cellCenter.x - 0.25f*m_dx, cellCenter.y + 0.25f*m_dx), // top left
					Vec2(cellCenter.x + 0.25f*m_dx, cellCenter.y + 0.25f*m_dx), // top right
					Vec2(cellCenter.x + 0.25f*m_dx, cellCenter.y - 0.25f*m_dx), // bottom right
					Vec2(cellCenter.x - 0.25f*m_dx, cellCenter.y - 0.25f*m_dx) // bottom left
				};
				// cycle through subgrid to place all particles
				for (int k = 0; k < particlesPerCell; k++) {
					// randomly jitter from subgrid center
					// give a random factor from [-0.24, 0.24] multiplied by dx
					float jitterX = ((float)((rand() % 49) - 24) / 100.0f) * m_dx;
					float jitterY = ((float)((rand() % 49) - 24) / 100.0f) * m_dx;
					Vec2 pos(subCenters[i % 4].x + jitterX, subCenters[i % 4].y + jitterY);
					Vec2 vel(0.0f, 0.0f);
					particleList->push_back(Particle2D(pos, vel));
				}
			}
		}
	}
}

/*
Updates the label grid based on particle positions. Grid cells containing particles
are fluid, solids remain the same, and all other are air. 
*/
void FluidSolver2D::labelGrid() {
	// first clear grid labels (mark everything air, but leave solids)
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (m_label[i][j] != SOLID) {
				m_label[i][j] = AIR;
			}
		}
	}

	// mark any cell containing a particle FLUID
	for (int i = 0; i < m_particles->size(); i++) {
		// get cell containing the particle
		int *cell = getGridCellIndex(m_particles->at(i).pos, m_dx);
		m_label[cell[0]][cell[1]] = FLUID;
	}
}

/*
Transfers particle values to the grid for cells labeled as a fluid. 
*/
void FluidSolver2D::particlesToGrid() {
	// init vel grids with unknown label value
	initGridValues<float>(m_u, m_gridWidth + 1, m_gridHeight, VEL_UNKNOWN);
	initGridValues<float>(m_v, m_gridWidth, m_gridHeight + 1, VEL_UNKNOWN);

	// For each component of velocity in each fluid grid cell
	// we calculate weighted average of particles around it defined
	// by a kernel function and set this as the vel in the cell. 

	// structures to accumulate grid numerator and denominator of weighted average before divide
	Mat2Dd uNum = initGrid2D<double>(m_gridWidth+1, m_gridHeight);
	Mat2Dd uDen = initGrid2D<double>(m_gridWidth+1, m_gridHeight);
	Mat2Dd vNum = initGrid2D<double>(m_gridWidth, m_gridHeight+1);
	Mat2Dd vDen = initGrid2D<double>(m_gridWidth, m_gridHeight+1);

	// clear accumulators
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			if (j < m_gridHeight) {
				uNum[i][j] = 0.0;
				uDen[i][j] = 0.0;
			}
			if (i < m_gridWidth) {
				vNum[i][j] = 0.0;
				vDen[i][j] = 0.0;
			}
		}
	}


	// loop over particles and accumulate num and den at each grid point
	for (int p = 0; p < m_particles->size(); p++) {
		Particle2D curParticle = m_particles->at(p);
		for (int i = 0; i < m_gridWidth + 1; i++) {
			for (int j = 0; j < m_gridHeight + 1; j++) {
					if (j < m_gridHeight) {
						double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i - 0.5f, j, m_dx)));
						uNum[i][j] += curParticle.vel.x * kernel;
						uDen[i][j] += kernel;
					}
					if (i < m_gridWidth) {
						double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i, j - 0.5f, m_dx)));
						vNum[i][j] += curParticle.vel.y * kernel;
						vDen[i][j] += kernel;
					}
			}
		}
	}

	if (DEBUG) printGrid2D<double>(m_gridWidth + 1, m_gridHeight, uNum);
	if (DEBUG) printGrid2D<double>(m_gridWidth + 1, m_gridHeight, uDen);
	if (DEBUG) printGrid2D<double>(m_gridWidth, m_gridHeight + 1, vNum);
	if (DEBUG) printGrid2D<double>(m_gridWidth, m_gridHeight + 1, vDen);

	// additional pass over grid to divide and update actual velocities
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
				if (j < m_gridHeight) {
					if (uDen[i][j] != 0.0) {
						m_u[i][j] = uNum[i][j] / uDen[i][j];
					}
				}
				if (i < m_gridWidth) {
					if (vDen[i][j] != 0.0) {
						m_v[i][j] = vNum[i][j] / vDen[i][j];
					}
				}
		}
	}

	deleteGrid2D<double>(m_gridWidth+1, m_gridHeight, uNum);
	deleteGrid2D<double>(m_gridWidth+1, m_gridHeight, uDen);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight+1, vNum);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight+1, vDen);
}

/*
Extrapolates the data in fluid cells of the given grid out using a breadth-first
search technique.
Args:
grid - the grid with data to extrapolate
x, y - the grid dimensions
depth - the number of cells away from fluid cells to extrapolate to.
*/
void FluidSolver2D::extrapolateGridFluidData(Mat2Df grid, int x, int y, int depth) {
	// initialize marker array 
	Mat2Di d = initGrid2D<int>(x, y);
	// set d to 0 for known values, max int for unknown
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			if (grid[i][j] != VEL_UNKNOWN) {
				d[i][j] = 0;
			} else {
				d[i][j] = INT_MAX;
			}
		}
	}

	// define neighbors
	int numNeighbors = 8;
	int neighbors[8][2] = {
		{-1, 1}, // top left
		{-1, 0}, // middle left
		{-1, -1}, // bottom left
		{0, 1}, // top middle
		{0, -1}, // bottom middle
		{1, 1}, // top right
		{1, 0}, // middle right
		{1, -1} // bottom right
	};

	// initialize first wavefront
	std::vector<Vec2> W;
	int dim[2] = { x, y };
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			// current value is not known
			if (d[i][j] != 0) {
				int ind[2] = { i, j };
				if (!checkNeighbors(d, dim, ind, neighbors, numNeighbors, 0).empty()) {
					// a neighbor is known
					d[i][j] = 1;
					W.push_back(Vec2(i, j));
				}
			}
		}
	}

	// list of all wavefronts, only want to go through the given depth
	std::vector<std::vector<Vec2>> wavefronts;
	wavefronts.push_back(W);
	int curWave = 0;
	while (curWave < depth) {
		// get wavefront
		std::vector<Vec2> curW = wavefronts.at(curWave);
		// initialize next wavefront
		std::vector<Vec2> nextW;
		// go through current wave and extrapolate values
		for (int i = 0; i < curW.size(); i++) {
			Vec2 ind = curW.at(i);
			// average neighbors
			float avg = 0.0f;
			int numUsed = 0;
			for (int i = 0; i < numNeighbors; i++) {
				int offsetX = neighbors[i][0];
				int offsetY = neighbors[i][1];
				int neighborX = ind.x + offsetX;
				int neighborY = ind.y + offsetY;

				// make sure valid indices
				if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1])) {
					// only want to add to average if neighbor d is less than current d
					if (d[neighborX][neighborY] < d[(int)ind.x][(int)ind.y]) {
						avg += grid[neighborX][neighborY];
						numUsed++;
					} else if (d[neighborX][neighborY] == INT_MAX) {
						d[neighborX][neighborY] = d[(int)ind.x][(int)ind.y] + 1;
						nextW.push_back(Vec2(neighborX, neighborY));
					}
				}
			}

			avg /= numUsed;
			// set current value to average of neighbors
			grid[(int)ind.x][(int)ind.y] = avg;
		}

		// push next wave to list
		wavefronts.push_back(nextW);
		curWave++;
	}

	// clean up
	deleteGrid2D<int>(x, y, d);
	
}

/*
Saves a copy of the current velocity grids to be used on the FLIP updated
*/
void FluidSolver2D::saveVelocityGrids() {
	// save u grid
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			m_uSaved[i][j] = m_u[i][j];
		}
	}

	// save v grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			m_vSaved[i][j] = m_v[i][j];
		}
	}
}

/*
Applies the force of gravity to velocity field on the grid
*/
void FluidSolver2D::applyBodyForces() {
	// traverse all grid cells and apply force to each velocity component
	// The new velocity is calculated using forward euler
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			if (j < m_gridHeight) {
				// make sure we know the velocity
				if (m_u[i][j] != VEL_UNKNOWN) {
					// update u component
					m_u[i][j] += m_dt*GRAVITY.x;
				}
			}
			if (i < m_gridWidth) {
				if (m_v[i][j] != VEL_UNKNOWN) {
					// update v component
					m_v[i][j] += m_dt*GRAVITY.y;
				}
			}
		}
	}
}

/*
Solves for pressure using the current velocity field.
*/
void FluidSolver2D::pressureSolve() {
	// initialize all grids to solve for pressure
	// using double for more accuracy
	Mat2Dd rhs = initGrid2D<double>(m_gridWidth, m_gridHeight);
	constructRHS(rhs);
	if (DEBUG) printGrid2D<double>(m_gridWidth, m_gridHeight, rhs);
	Mat2Dd Adiag = initGrid2D<double>(m_gridWidth, m_gridHeight);
	Mat2Dd Ax = initGrid2D<double>(m_gridWidth, m_gridHeight);
	Mat2Dd Ay = initGrid2D<double>(m_gridWidth, m_gridHeight);
	constructA(Adiag, Ax, Ay);
	if (DEBUG) printGrid2D<double>(m_gridWidth, m_gridHeight, Adiag);
	if (DEBUG) printGrid2D<double>(m_gridWidth, m_gridHeight, Ax);
	if (DEBUG) printGrid2D<double>(m_gridWidth, m_gridHeight, Ay);
	Mat2Dd precon = initGrid2D<double>(m_gridWidth, m_gridHeight);
	constructPrecon(precon, Adiag, Ax, Ay);
	if (DEBUG) printGrid2D<double>(m_gridWidth, m_gridHeight, precon);

	// solve for pressure using PCG
	PCG(Adiag, Ax, Ay, rhs, precon);

	// cleanup
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, rhs);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, Adiag);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, Ax);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, Ay);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, precon);
}

/*
Applies the pressure force to the current velocity field.
*/
void FluidSolver2D::applyPressure() {
	float scale = m_dt / (FLUID_DENSITY * m_dx);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			// update u
			if (i - 1 >= 0) {
				if (m_label[i - 1][j] == FLUID || m_label[i][j] == FLUID) {
					if (m_label[i - 1][j] == SOLID || m_label[i][j] == SOLID) {
						// TODO add solid velocities
						m_u[i][j] = 0.0f; // usolid[i][j]
					} else {
						m_u[i][j] -= scale * (m_p[i][j] - m_p[i - 1][j]);
					}
				} else {
					m_u[i][j] = VEL_UNKNOWN;
				}
			} else {
				// edge of grid, keep the same velocity
			}

			// update v
			if (j - 1 >= 0) {
				if (m_label[i][j - 1] == FLUID || m_label[i][j] == FLUID) {
					if (m_label[i][j - 1] == SOLID || m_label[i][j] == SOLID) {
						// TODO add solid velocities
						m_v[i][j] = 0.0f; // vsolid[i][j]
					}
					else {
						m_v[i][j] -= scale * (m_p[i][j] - m_p[i][j-1]);
					}
				} else {
					m_v[i][j] = VEL_UNKNOWN;
				}
			} else {
				// edge of grid, keep the same velocity
			}
		}
	}
}

/*
Transfer the velocities from the grid back to the particles. This is done
with a PIC/FLIP mix, where the PIC update has a weight of the given alpha.
Args:
alpha - weight in the update for PIC, should be in [0, 1]. Then the FLIP update is weighted (1 - alpha).
*/
void FluidSolver2D::gridToParticles(float alpha) {
	// build grid for change in velocity to use for FLIP update
	Mat2Df duGrid = initGrid2D<float>(m_gridWidth + 1, m_gridHeight);
	Mat2Df dvGrid = initGrid2D<float>(m_gridWidth, m_gridHeight + 1);
	// calc u grid
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			duGrid[i][j] = m_u[i][j] - m_uSaved[i][j];
		}
	}
	// calc v grid
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			dvGrid[i][j] = m_v[i][j] - m_vSaved[i][j];
		}
	}

	// go through particles and interpolate each velocity component
	// the update is a PIC/FLIP mix weighted with alpha
	// alpha = 1.0 is entirely PIC, alpha = 0.0 is all FLIP
	for (int i = 0; i < m_particles->size(); i++) {
		Particle2D *curParticle = &(m_particles->at(i));
		Vec2 picInterp = interpVel(m_u, m_v, curParticle->pos);
		Vec2 flipInterp = interpVel(duGrid, dvGrid, curParticle->pos);
		// u_new = alpha * interp(u_gridNew, x_p) + (1 - alpha) * (u_pOld + interp(u_dGrid, x_p))
		curParticle->vel = add(scale(picInterp, alpha), scale(add(curParticle->vel, flipInterp), 1.0f - alpha));
	}

	deleteGrid2D<float>(m_gridWidth + 1, m_gridHeight, duGrid);
	deleteGrid2D<float>(m_gridWidth, m_gridHeight + 1, dvGrid);
}

/*
Advects the particles through the current velocity field using a Runge-Kutta 3 method.
This uses substeps so that the particles are never advected greater than C*dx in a single
substep.
Args:
C - the maximum number of grid cells a particle should move when advected. This helps define substep sizes. 
*/
void FluidSolver2D::advectParticles(int C) {
	for (int i = 0; i < m_particles->size(); i++) {
		Particle2D *curParticle = &(m_particles->at(i));
		float subTime = 0;
		bool finished = false;
		//float dT = m_dt / 4.999f;
		while (!finished) {
			Vec2 curVel = interpVel(m_u, m_v, curParticle->pos);

			// calc max substep size
			float dT = (C * m_dx) / (norm(curVel) + FLT_MIN);
			// update substep time so we don't go past normal time step
			if (subTime + dT >= m_dt) {
				dT = m_dt - subTime;
				finished = true;
			} else if (subTime + 2 * dT >= m_dt) {
				dT = 0.5f * (m_dt - subTime);
			}

			RK3(curParticle, curVel, dT, m_u, m_v);
			subTime += dT;

			if (curParticle->pos.x < 0 || curParticle->pos.y < 0 || isnan(curParticle->pos.x) || isnan(curParticle->pos.y)) {
				// there's been an error in RK3, just skip it
				std::cout << "RK3 error...skipping particle" << std::endl;
				break;
			}

			int *cell = getGridCellIndex(curParticle->pos, m_dx);
			int j = cell[0];
			int k = cell[1];
			if (m_label[j][k] == SOLID) {
				//std::cout << "Advected into SOLID, projecting back!\n";
				if (!projectParticle(curParticle, m_dx / 4.0f)) {
					std::cout << "RK3 error...skipping particle" << std::endl;
					break;
				}
			}

		}
	}
}

/*
Finds particles that have traveled into solid boundaries and projects them back into the fluid region. 
They will be put at the closest boundary + the given dx into the fluid.
Args
dx - the amount to project stray particles away from the wall.
*/
void FluidSolver2D::cleanupParticles(float dx) {
	int i = 0;
	bool finished = false;
	int numDeleted = 0;
	while(!finished && m_particles->size() > 0) {
		int *cell = getGridCellIndex(m_particles->at(i).pos, m_dx);
		int ind[2] = { cell[0], cell[1] };
		// if either of cells are negative or greater than sim dimensions it has left sim area
		if (ind[0] < 0 || ind[1] < 0 || ind[0] >= m_gridWidth || ind[1] >= m_gridHeight || isnan(m_particles->at(i).pos.x) || isnan(m_particles->at(i).pos.y)) {
			m_particles->erase(m_particles->begin() + i);
			numDeleted++;
			if (i >= m_particles->size()) {
				finished = true;
			}
		} else if (m_label[ind[0]][ind[1]] == SOLID) {
			// project back into fluid
			bool success = projectParticle(&(m_particles->at(i)), dx);
			if (!success) {
				// no near fluid, just delete
				m_particles->erase(m_particles->begin() + i);
				numDeleted++;
				if (i >= m_particles->size()) {
					finished = true;
				}
			} 
		} else {
			i++;
			if (i >= m_particles->size()) {
				finished = true;
			}
		}
	}

	if (DEBUG) std::cout << "Removed " << numDeleted << " particles from sim.\n";
}



//----------------------------------------------------------------------
// Private Helper Functions
//----------------------------------------------------------------------

/*
Populates the given grid with the value given.
Args:
grid - the grid to fill
x/y - the dimensions of the grid
value - the value to fill it with
*/
template <typename T>
void FluidSolver2D::initGridValues(T** grid, int x, int y, T value) {
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			grid[i][j] = value;
		}
	}
}

/*
Checks neighbors of the given index in the given grid for the given value. Returns a vector
of neighbor indices (row index in the given neighbor array) that have the given value.
Args
grid - the 2D grid to look in
dim - the grid dimensions [x y]
index - the (i, j) index of the cell to look around
neighbors - the definintion of neighbors, this is a n x 2 array where each row is a pair of offsets from the index of interest
numNeighbors - the number of neighbors
value - the value to look for
*/
std::vector<int> FluidSolver2D::checkNeighbors(Mat2Di grid, int dim[2], int index[2], int neighbors[][2], int numNeighbors, int value) {
	std::vector<int> neighborsTrue;
	for (int i = 0; i < numNeighbors; i++) {
		int offsetX = neighbors[i][0];
		int offsetY = neighbors[i][1];
		int neighborX = index[0] + offsetX;
		int neighborY = index[1] + offsetY;

		// make sure valid indices
		if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1])) {
			if (grid[neighborX][neighborY] == value) {
				neighborsTrue.push_back(i);
			}
		}
	}

	return neighborsTrue;
}


/*
Returns the value of the trilinear hat function for the given
distance (x, y).
*/
double FluidSolver2D::trilinearHatKernel(SimUtil::Vec2 dist) {
	return hatFunction(dist.x / m_dx) * hatFunction(dist.y / m_dx);
}

/*
Calculates the value of hat function for the given r. 
*/
double FluidSolver2D::hatFunction(double r) {
	double rAbs = abs(r);
	if (rAbs <= 1) {
		return 1.0 - rAbs;
	} else {
		return 0.0;
	}
}


/*
Returns the value of the quadratic B-spline function for the
given distance (x, y).
*/
double FluidSolver2D::quadBSplineKernel(SimUtil::Vec2) {
	// TODO
	return 0;
}

/*
Interpolates the value in the given velocity grid at the given position using bilinear interpolation.
Returns velocity unkown if position is not on simulation grid. 
Args:
uGrid - the u component grid to interpolate from
vGrid - the v component grid to interpolate from
pos - the position to interpolate at
*/
Vec2 FluidSolver2D::interpVel(SimUtil::Mat2Df uGrid, SimUtil::Mat2Df vGrid, Vec2 pos) {
	// get grid cell containing position
	int *cell = getGridCellIndex(pos, m_dx);
	int i = cell[0];
	int j = cell[1];
	// make sure this is a valid index
	if (i >= 0 && i < m_gridWidth && j >= 0 && j < m_gridHeight) {
		// get positions of u and v component stored on each side of cell
		Vec2 cellLoc = getGridCellPosition(i, j, m_dx);
		float offset = m_dx / 2.0f;
		float x1 = cellLoc.x - offset;
		float x2 = cellLoc.x + offset;
		float y1 = cellLoc.y - offset;
		float y2 = cellLoc.y + offset;
		// get actual values at these positions
		float u1 = uGrid[i][j];
		float u2 = uGrid[i + 1][j];
		float v1 = vGrid[i][j];
		float v2 = vGrid[i][j + 1];

		// the interpolated values
		float u = ((x2 - pos.x) / (x2 - x1)) * u1 + ((pos.x - x1) / (x2 - x1)) * u2;
		float v = ((y2 - pos.y) / (y2 - y1)) * v1 + ((pos.y - y1) / (y2 - y1)) * v2;
		return Vec2(u, v);
	} else {
		return Vec2(VEL_UNKNOWN, VEL_UNKNOWN);
	}
}

/*
Advects a particle using Runge-Kutta 3 through the given velocity field.
Args:
particle - the particle to advect
initVel - the particles initial velocity in the current field, can leave UNKNOWN
dt - the time step
uGrid/vGrid - the velocity grids to advect through
*/
void FluidSolver2D::RK3(SimUtil::Particle2D *particle, SimUtil::Vec2 initVel, float dt, SimUtil::Mat2Df uGrid, SimUtil::Mat2Df vGrid) {
	if (initVel.x == VEL_UNKNOWN && initVel.y == VEL_UNKNOWN) {
		initVel = interpVel(uGrid, vGrid, particle->pos);
	}

	Vec2 k1 = initVel;
	Vec2 k2 = interpVel(uGrid, vGrid, add(particle->pos, scale(k1, 0.5f*dt)));
	Vec2 k3 = interpVel(uGrid, vGrid, add(particle->pos, scale(k2, 0.75f*dt)));
	k1 = scale(k1, (2.0f / 9.0f)*dt);
	k2 = scale(k2, (3.0f / 9.0f)*dt);
	k3 = scale(k3, (4.0f / 9.0f)*dt);

	particle->pos = add(particle->pos, add(k1, add(k2, k3)));
}

/*
Sets up the right hand side of the system to solve for pressure. This is the negative
divergence at each cell center modified to account for the velocity of solids at boundaries.
Args:
rhs - the grid to use for the RHS
*/
void FluidSolver2D::constructRHS(Mat2Dd rhs) {
	// initialize to 0
	initGridValues<double>(rhs, m_gridWidth, m_gridHeight, 0.0);
	// calculate negative divergence
	double scale = 1.0f / m_dx;
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (isFluid(i, j)) {
				rhs[i][j] = -scale * (m_u[i + 1][j] - m_u[i][j] + m_v[i][j + 1] - m_v[i][j]);
				// if it's on boundary must update to consider solid velocity
				// TODO create actual solid velocity grids, for right now just 0
				if (m_label[i - 1][j] == SOLID) {
					rhs[i][j] -= scale * (m_u[i][j] - 0.0f); //m_usolid[i][j]
				}
				if (m_label[i + 1][j] == SOLID) {
					rhs[i][j] += scale * (m_u[i+1][j] - 0.0f); //m_usolid[i+1][j]
				}
				if (m_label[i][j - 1] == SOLID) {
					rhs[i][j] -= scale * (m_v[i][j] - 0.0f); //m_vsolid[i][j]
				}
				if (m_label[i][j + 1] == SOLID) {
					rhs[i][j] += scale * (m_v[i][j+1] - 0.0f); //m_vsolid[i][j+1]
				}
			}
		}
	}
}

/*
Constructs the A matrix for the system to solve for pressure. This a sparse coefficient matrix
for the pressure terms, stored in 3 separate grids. If index i, j, k is not a fluid cell, then
it is 0.0 in all 3 grids that store the matix. 
Args:
Adiag - grid to store the diagonal of the matrix in.
Ax - grid to store the coefficients for pressure in the (i+1) cell for each grid cell with x index i
Ay - grid to store the coefficients for pressure in the (j+1) cell for each grid cell with y index j
*/
void FluidSolver2D::constructA(Mat2Dd Adiag, Mat2Dd Ax, Mat2Dd Ay) {
	// clear to all zeros so can increment
	initGridValues<double>(Adiag, m_gridWidth, m_gridHeight, 0.0);
	initGridValues<double>(Ax, m_gridWidth, m_gridHeight, 0.0);
	initGridValues<double>(Ay, m_gridWidth, m_gridHeight, 0.0);

	// populate with coefficients for pressure unknowns
	double scale = m_dt / (FLUID_DENSITY * m_dx * m_dx);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (isFluid(i, j)) {
				// handle negative x neighbor
				if (m_label[i - 1][j] == FLUID || m_label[i - 1][j] == AIR) {
					Adiag[i][j] += scale;
				}
				// handle positive x neighbor
				if (m_label[i + 1][j] == FLUID) {
					Adiag[i][j] += scale;
					Ax[i][j] = -scale;
				} else if (m_label[i + 1][j] == AIR) {
					Adiag[i][j] += scale;
				}
				// handle negative y neighbor
				if (m_label[i][j - 1] == FLUID || m_label[i][j - 1] == AIR) {
					Adiag[i][j] += scale;
				}
				// handle positive y neighbor
				if (m_label[i][j + 1] == FLUID) {
					Adiag[i][j] += scale;
					Ay[i][j] = -scale;
				} else if (m_label[i][j + 1] == AIR) {
					Adiag[i][j] += scale;
				}
			}
		}
	}
}

/*
Constructs the preconditioner used when performing the preconditioned conjugate gradient (PCG)
algorithm to solve for pressure.
Args:
precon - grid to store the preconditioner in
Adiag, Ax, Ay - the grids that make up the A coefficient matrix
*/
void FluidSolver2D::constructPrecon(Mat2Dd precon, Mat2Dd Adiag, Mat2Dd Ax, Mat2Dd Ay) {
	initGridValues<double>(precon, m_gridWidth, m_gridHeight, 0.0);

	// tuning constant
	double tau = 0.97;
	// safety constant
	double sigma = 0.25;

	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (isFluid(i, j)) {
				double Adiag_ij = Adiag[i][j];
				double Ax_im1j = 0.0;
				double Ax_ijm1 = 0.0;
				double Ay_ijm1 = 0.0;
				double Ay_im1j = 0.0;
				double precon_im1j = 0.0;
				double precon_ijm1 = 0.0;
				// want to stay at zero if off the grid
				// non-fluid entries in A are already 0
				if (i - 1 >= 0 && i - 1 < m_gridWidth) {
					if (isFluid(i - 1, j)) {
						Ax_im1j = Ax[i - 1][j];
						Ay_im1j = Ay[i - 1][j];
						precon_im1j = precon[i - 1][j];
					}
				}
				if (j - 1 >= 0 && j - 1 < m_gridHeight) {
					if (isFluid(i, j - 1)) {
						Ax_ijm1 = Ax[i][j - 1];
						Ay_ijm1 = Ay[i][j - 1];
						precon_ijm1 = precon[i][j - 1];
					}
				}

				double e = Adiag_ij - pow(Ax_im1j * precon_im1j, 2.0)
							- pow(Ay_ijm1 * precon_ijm1, 2.0)
							- tau * (
								Ax_im1j * Ay_im1j * pow(precon_im1j, 2.0)
								+ Ay_ijm1 * Ax_ijm1 * pow(precon_ijm1, 2.0)
								);

				if (e < (sigma * Adiag_ij)) {
					e = Adiag_ij;
				}

				precon[i][j] = 1.0 / sqrt(e);
			}
		}
	}
}

/*
Performs the Modified Incomplete Cholesky Conjugate Gradient Level Zero 
(preconditioned conjugate gradient, PCG) algorithm to solve the linear system
Ap = b for p. The result are placed on the pressure grid.
Args:
Adiag, Ax, Ay - the grids that make up the A coefficient matrix
b - the right hand side of the equation
precon - the preconditioner to use
*/
void FluidSolver2D::PCG(Mat2Dd Adiag, Mat2Dd Ax, Mat2Dd Ay, Mat2Dd b, Mat2Dd precon) {
	// reset pressure vec to 0
	initGridValues<float>(m_p, m_gridWidth, m_gridHeight, 0.0f);
	// residiual vector is b to start
	Mat2Dd r = initGrid2D<double>(m_gridWidth, m_gridHeight);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			r[i][j] = b[i][j];
		}
	}

	// check if r = 0
	bool r0 = true;
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (r[i][j] != 0) {
				r0 = false;
				break;
			}
		}
	}

	if (r0) {
		std::cout << "Did not run PCG b/c of 0 residual to start.\n";
		deleteGrid2D<double>(m_gridWidth, m_gridHeight, r);
		return;
	}

	// auxiliary vector
	Mat2Dd z = initGrid2D<double>(m_gridWidth, m_gridHeight);
	// search vector
	Mat2Dd s = initGrid2D<double>(m_gridWidth, m_gridHeight);

	// initialize auxiliary and search
	applyPrecon(z, r, precon, Adiag, Ax, Ay);
	// search is set to auxiliary to start
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			s[i][j] = z[i][j];
		}
	}

	double sigma = dot(z, r, m_gridWidth, m_gridHeight);

	// start the main iteration until tolerance is reach or max out iterations
	bool converged = false;
	for (int iters = 0; iters < PCG_MAX_ITERS; iters++) {
		applyA(z, s, Adiag, Ax, Ay);
		double alpha = sigma / dot(z, s, m_gridWidth, m_gridHeight);
		// update pressure and residual
		for (int i = 0; i < m_gridWidth; i++) {
			for (int j = 0; j < m_gridHeight; j++) {
				m_p[i][j] += (alpha * s[i][j]);
				r[i][j] -= (alpha * z[i][j]);
			}
		}
		// check if we're under the tolerance
		if (max(r, m_gridWidth, m_gridHeight) <= PCG_TOL) {
			if (DEBUG) std::cout << "PCG converged after " << iters << " iterations.\n";
			converged = true;
			break;
		}
		// otherwise new auxiliary vector
		applyPrecon(z, r, precon, Adiag, Ax, Ay);
		double newSigma = dot(z, r, m_gridWidth, m_gridHeight);
		double beta = newSigma / sigma;
		// update search vector
		for (int i = 0; i < m_gridWidth; i++) {
			for (int j = 0; j < m_gridHeight; j++) {
				s[i][j] = z[i][j] + (beta * s[i][j]);
			}
		}
		// update sigma
		sigma = newSigma;
	}

	if (!converged) {
		std::cout << "PCG did not converge, stopped after " << PCG_MAX_ITERS << " iterations!\n";
	}

	// cleanup 
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, r);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, z);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, s);
}

/*
Applies the given preconditioner to the given vector and places the result in given "vector" z.
Args:
z - the vector (2D grid) to place the result in
r - the vector (2D grid) to multiply the preconditioner by
precon - the preconditioner
Adiag, Ax, Ay - the grids that make up the coefficient matrix
*/
void FluidSolver2D::applyPrecon(Mat2Dd z, Mat2Dd r, Mat2Dd precon, Mat2Dd Adiag, Mat2Dd Ax, Mat2Dd Ay) {
	// first solve Lq = r
	Mat2Dd q = initGrid2D<double>(m_gridWidth, m_gridHeight);
	initGridValues<double>(q, m_gridWidth, m_gridHeight, 0.0);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (isFluid(i, j)) {
				double Ax_im1j = 0.0;
				double Ay_ijm1 = 0.0;
				double precon_im1j = 0.0;
				double precon_ijm1 = 0.0;
				double q_im1j = 0.0;
				double q_ijm1 = 0.0;

				if (i - 1 >= 0 && i - 1 < m_gridWidth) {
					if (isFluid(i - 1, j)) {
						Ax_im1j = Ax[i - 1][j];
						precon_im1j = precon[i - 1][j];
						q_im1j = q[i - 1][j];
					}
				}
				if (j - 1 >= 0 && j - 1 < m_gridHeight) {
					if (isFluid(i, j - 1)) {
						Ay_ijm1 = Ay[i][j - 1];
						precon_ijm1 = precon[i][j - 1];
						q_ijm1 = q[i][j - 1];
					}
				}

				double t = r[i][j] - (Ax_im1j * precon_im1j * q_im1j)
								   - (Ay_ijm1 * precon_ijm1 * q_ijm1);

				q[i][j] = t * precon[i][j];
			}
		}
	}

	// now solve L^T z = q
	initGridValues<double>(z, m_gridWidth, m_gridHeight, 0.0);
	for (int i = m_gridWidth - 1; i >= 0; i--) {
		for (int j = m_gridHeight; j >= 0; j--) {
			if (isFluid(i, j)) {
				double Ax_ij = Ax[i][j];
				double Ay_ij = Ay[i][j];
				double precon_ij = precon[i][j];
				double z_ip1j = 0.0;
				double z_ijp1 = 0.0;

				if (i + 1 >= 0 && i + 1 < m_gridWidth) {
					if (isFluid(i + 1, j)) {
						z_ip1j = z[i + 1][j];
					}
				}
				if (j + 1 >= 0 && j + 1 < m_gridHeight) {
					if (isFluid(i, j + 1)) {
						z_ijp1 = z[i][j + 1];
					}
				}

				double t = q[i][j] - (Ax_ij * precon_ij * z_ip1j)
								   - (Ay_ij * precon_ij * z_ijp1);

				z[i][j] = t * precon_ij;
			}
		}
	}

	// cleanup
	deleteGrid2D<double>(m_gridWidth, m_gridHeight, q);
}

/*
Multiplies A by the given vector (a 2D grid).
Args:
z - the vector (2D grid) to place the result in
s - the vector (2D grid) to multiply A by
Adiag, Ax, Ay - the grids that make up the coefficient matrix A
*/
void FluidSolver2D::applyA(Mat2Dd z, Mat2Dd s, Mat2Dd Adiag, Mat2Dd Ax, Mat2Dd Ay) {
	initGridValues<double>(z, m_gridWidth, m_gridHeight, 0.0);
	for (int i = 0; i < m_gridWidth; i++) {
		for (int j = 0; j < m_gridHeight; j++) {
			if (isFluid(i, j)) {
				z[i][j] = Adiag[i][j] * s[i][j]
					+ Ax[i][j] * s[i + 1][j]
					+ Ay[i][j] * s[i][j + 1];
				if (i - 1 >= 0 && i - 1 < m_gridWidth) {
					z[i][j] += Ax[i - 1][j] * s[i - 1][j];
				}
				if (j - 1 >= 0 && j - 1 < m_gridHeight) {
					z[i][j] += Ay[i][j - 1] * s[i][j - 1];
				}
			}
		}
	}
}

/*
Projects a particle within a solid back into the closest fluid or air. Returns true if
successful and false otherwise.
Args
particle - the particle to project. 
dx - the amount to project the particle from the solid wall
*/
bool FluidSolver2D::projectParticle(Particle2D *particle, float dx) {
	// project back into fluid
	// find neighbors that are fluid
	// define neighbors
	int numNeighbors = 8;
	int neighbors[8][2] = {
		{ -1, 1 }, // top left
		{ -1, 0 }, // middle left
		{ -1, -1 }, // bottom left
		{ 0, 1 }, // top middle
		{ 0, -1 }, // bottom middle
		{ 1, 1 }, // top right
		{ 1, 0 }, // middle right
		{ 1, -1 } // bottom right
	};
	int dim[2] = { m_gridWidth, m_gridHeight };
	int *cell = getGridCellIndex(particle->pos, m_dx);
	int index[2] = { cell[0], cell[1] };
	// get neighbors that are fluid
	std::vector<int> neighborInd = checkNeighbors(m_label, dim, index, neighbors, numNeighbors, FLUID);
	if (neighborInd.size() == 0) {
		// try with air
		neighborInd = checkNeighbors(m_label, dim, index, neighbors, numNeighbors, AIR);
	}
	// find closest to particle
	int closestInd = -1;
	float closestDist = std::numeric_limits<float>::max();
	Vec2 closestVec(0.0f, 0.0f);
	for (int j = 0; j < neighborInd.size(); j++) {
		// get vec from particle to neighbor ind
		int ind[2] = { cell[0] + neighbors[neighborInd.at(j)][0], cell[1] + neighbors[neighborInd.at(j)][1] };
		Vec2 cellPos = getGridCellPosition(ind[0], ind[1], m_dx);
		Vec2 distVec = sub(cellPos, particle->pos);
		float dist = norm(distVec);
		if (dist < closestDist) {
			closestDist = dist;
			closestInd = neighborInd.at(j);
			closestVec = distVec;
		}
	}

	if (closestInd == -1) {
		return false;
	}
	else {
		// project different ways based on where closest neighbor is
		// also make sure to only project the amount given
		Vec2 projectVec(0.0f, 0.0f);
		if (closestInd == 1) { // middle left
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
		}
		else if (closestInd == 3) { // top middle
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 4) { // bottom middle
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
		}
		else if (closestInd == 6) { // middle right
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 5) { // top right
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 0) { // top left
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (dx - (m_dx / 2.0f));
		}
		else if (closestInd == 2) { // bottom left
			projectVec.x = closestVec.x + (-dx + (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
		}
		else if (closestInd == 7) { // bottom right
			projectVec.x = closestVec.x + (dx - (m_dx / 2.0f));
			projectVec.y = closestVec.y + (-dx + (m_dx / 2.0f));
		}

		particle->pos = add(particle->pos, projectVec);

		return true;
	}
}


/*
Determines if the given grid cell is considered a fluid based on the label grid. Also takes
into account velocity components on the edge of the grid. For example, if the index passed in
is one cell outside the label grid, it is assumed to be a velocity component index, and whether the cell is
fluid or not is determined by the cell that it borders. Otherwise false is returned.
Args
i - x cell index
j - y cell index
*/
bool FluidSolver2D::isFluid(int i, int j) {
	bool isFluid = false;
	// see if velocity on edge of grid
	// if it is we can't check the label at that point, must check the one previous
	if (i == m_gridWidth || j == m_gridHeight) {
		// i and j should never both be out of label range
		// should only happen in one dimension because of vel comp grids
		if (i == m_gridWidth && j == m_gridHeight) {
			isFluid = false;
		}
		else if (i == m_gridWidth) {
			if (m_label[i - 1][j] == FLUID) {
				isFluid = true;
			}
		}
		else if (j == m_gridHeight) {
			if (m_label[i][j - 1] == FLUID) {
				isFluid = true;
			}
		}
	}
	else if (m_label[i][j] == FLUID) {
		isFluid = true;
	}

	return isFluid;
}

