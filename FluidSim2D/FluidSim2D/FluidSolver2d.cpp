#include <iostream>
#include <exception>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "FluidSolver2d.h"
#include "SimUtil.h"

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
	deleteGrid2D<float>(m_gridWidth, m_gridHeight + 1, m_v);

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
	m_v = initGrid2D<float>(m_gridWidth, m_gridHeight + 1);

	// read in initial geometry to populate label grid
	readInGeom(m_gridWidth, m_gridHeight, initialGeometryFile, m_label);
	//printMat2D<int>(m_gridHeight, m_gridWidth, m_label);

	// seed particles using label grid
	seedParticles(PARTICLES_PER_CELL, m_particles);
}

void FluidSolver2D::step() {
	// update the grid labels
	labelGrid();
	// transfer particle vel to grid
	particlesToGrid();
	// extrapolate fluid data out one cell for accurate divergence calculations
	extrapolateGridFluidData(1);

}

void FluidSolver2D::saveParticleData(std::ofstream *particleOut) {
	if (particleOut->is_open()) {
		// print out all particle data on same line, each pos separated by ","
		size_t numParticles = m_particles->size();
		for (int i = 0; i < numParticles - 1; i++) {
			(*particleOut) << m_particles->at(i).pos.x << " " << m_particles->at(i).pos.y << ",";
		}
		(*particleOut) << m_particles->at(numParticles - 1).pos.x << " " << m_particles->at(numParticles - 1).pos.y << "\n";
	}
}



//----------------------------------------------------------------------
// Private Functions
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
					// give a random factor from [-0.25, 0.25] multiplied by dx
					float jitterX = ((float)((rand() % 51) - 25) / 100.0f) * m_dx;
					float jitterY = ((float)((rand() % 51) - 25) / 100.0f) * m_dx;
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
				// see if velocity on edge of grid
				// if it is we can't check the label at that point, must check the one previous
				bool isFluid = false;
				if (i == m_gridWidth || j == m_gridHeight) {
					if (i == m_gridWidth && j == m_gridHeight) {
						if (m_label[i - 1][j - 1] == FLUID) {
							isFluid = true;
						}
					} else if (i == m_gridWidth) {
						if (m_label[i - 1][j] == FLUID) {
							isFluid = true;
						}
					} else if (j == m_gridHeight) {
						if (m_label[i][j - 1] == FLUID) {
							isFluid = true;
						}
					}
				} else if (m_label[i][j] == FLUID) {
					isFluid = true;
				}

				if (isFluid) {
					double kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i - 0.5f, j, m_dx)));
					uNum[i][j] += curParticle.vel.x * kernel;
					uDen[i][j] += kernel;

					kernel = trilinearHatKernel(sub(curParticle.pos, getGridCellPosition(i, j - 0.5f, m_dx)));
					vNum[i][j] += curParticle.vel.y * kernel;
					vDen[i][j] += kernel;
				}
			}
		}
	}

	// additional pass over grid to divide and update actual velocities
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			bool isFluid = false;
			if (i == m_gridWidth || j == m_gridHeight) {
				if (i == m_gridWidth && j == m_gridHeight) {
					if (m_label[i - 1][j - 1] == FLUID) {
						isFluid = true;
					}
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

			if (isFluid) {
				m_u[i][j] = uNum[i][j] / uDen[i][j];
				m_v[i][j] = vNum[i][j] / vDen[i][j];
			}
		}
	}

	deleteGrid2D<double>(m_gridWidth+1, m_gridHeight, uNum);
	deleteGrid2D<double>(m_gridWidth+1, m_gridHeight, uDen);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight+1, vNum);
	deleteGrid2D<double>(m_gridWidth, m_gridHeight+1, vDen);
}

/*
Extrapolates the data (velocity) in fluid cells out using a breadth-first
search technique.
Args:
depth - the number of cells away from fluid cells to extrapolate to.
*/
void FluidSolver2D::extrapolateGridFluidData(int depth) {
	// TODO
	return;
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
	if (r <= 1) {
		return 1.0 - r;
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



