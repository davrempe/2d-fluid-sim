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
	deleteMat2D<int>(m_gridHeight, m_gridWidth, m_label);
	deleteMat2D<float>(m_gridHeight, m_gridWidth, m_p);
	deleteMat2D<float>(m_gridHeight, m_gridWidth + 1, m_u);
	deleteMat2D<float>(m_gridHeight + 1, m_gridWidth, m_v);

	delete m_particles;
}

//----------------------------------------------------------------------
// Public Functions
//----------------------------------------------------------------------

void FluidSolver2D::init(std::string initialGeometryFile){
	// set up the grid for simulation by initializing arrays
	m_label = initMat2D<int>(m_gridHeight, m_gridWidth);
	m_p = initMat2D<float>(m_gridHeight, m_gridWidth);
	m_u = initMat2D<float>(m_gridHeight, m_gridWidth + 1);
	m_v = initMat2D<float>(m_gridHeight + 1, m_gridWidth);

	// read in initial geometry to populate label grid
	readInGeom(m_gridWidth, m_gridHeight, initialGeometryFile, m_label);
	//printMat2D<int>(m_gridHeight, m_gridWidth, m_label);

	// seed particles using label grid
	seedParticles(PARTICLES_PER_CELL, m_particles);
}

void FluidSolver2D::step() {

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
	for (int i = 0; i < m_gridHeight; i++) {
		for (int j = 0; j < m_gridWidth; j++) {
			if (m_label[i][j] == SimUtil::FLUID) {
				// seed randomly in 2x2 subgrid of the cell
				Vec2 cellCenter = getCellLocation(i, j, m_dx);
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



