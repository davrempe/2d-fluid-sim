#include <iostream>
#include <exception>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <climits>

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
	extrapolateGridFluidData(m_u, m_gridWidth + 1, m_gridHeight, 1);
	extrapolateGridFluidData(m_v, m_gridWidth, m_gridHeight + 1, 1);
	
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
				if (isFluid(i, j)) {
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
	}

	// additional pass over grid to divide and update actual velocities
	for (int i = 0; i < m_gridWidth + 1; i++) {
		for (int j = 0; j < m_gridHeight + 1; j++) {
			if (isFluid(i, j)) {
				if (j < m_gridHeight) {
					m_u[i][j] = uNum[i][j] / uDen[i][j];
				}
				if (i < m_gridWidth) {
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
			if (isFluid(i, j)) {
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
				if (checkNeighbors(d, dim, ind, neighbors, numNeighbors, 0)) {
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
	
}

/*
Checks neighbors of the given index in the given grid for the given value. Returns true if
any neighbor has that value, and false otherwise.
Args
grid - the 2D grid to look in
dim - the grid dimensions [x y]
index - the (i, j) index of the cell to look around
neighbors - the definintion of neighbors, this is a n x 2 array where each row is a pair of offsets from the index of interest
numNeighbors - the number of neighbors
value - the value to look for
*/
bool FluidSolver2D::checkNeighbors(Mat2Di grid, int dim[2], int index[2], int neighbors[][2], int numNeighbors, int value) {
	for (int i = 0; i < numNeighbors; i++) {
		int offsetX = neighbors[i][0];
		int offsetY = neighbors[i][1];
		int neighborX = index[0] + offsetX;
		int neighborY = index[1] + offsetY;

		// make sure valid indices
		if ((neighborX >= 0 && neighborX < dim[0]) && (neighborY >= 0 && neighborY < dim[1])) {
			if (grid[neighborX][neighborY] == value) {
				return true;
			}
		}
	}

	return false;
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

