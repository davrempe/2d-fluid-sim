#include <fstream>
#include <iostream>
#include <exception>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "FluidSolver2d.h"
#include "SimUtil.h"


FluidSolver2D::FluidSolver2D(int width, int height, float dx, float dt){
	m_gridWidth = width;
	m_gridHeight = height;
	m_dx = dx;
	m_dt = dt;
}

FluidSolver2D::~FluidSolver2D(){
	// clean up grid
	SimUtil::deleteMat2D<int>(m_gridHeight, m_gridWidth, m_label);
	SimUtil::deleteMat2D<float>(m_gridHeight, m_gridWidth, m_p);
	SimUtil::deleteMat2D<float>(m_gridHeight, m_gridWidth + 1, m_u);
	SimUtil::deleteMat2D<float>(m_gridHeight + 1, m_gridWidth, m_v);

	delete m_particles;
}


void FluidSolver2D::init(std::string initialGeometryFile){
	// set up the grid for simulation by initializing arrays
	m_label = SimUtil::initMat2D<int>(m_gridHeight, m_gridWidth);
	m_p = SimUtil::initMat2D<float>(m_gridHeight, m_gridWidth);
	m_u = SimUtil::initMat2D<float>(m_gridHeight, m_gridWidth + 1);
	m_v = SimUtil::initMat2D<float>(m_gridHeight + 1, m_gridWidth);

	// read in initial geometry to populate label grid
	readInGeom(m_gridWidth, m_gridHeight, initialGeometryFile, m_label);
	//SimUtil::printMat2D<int>(m_gridHeight, m_gridWidth, m_label);

	// seed particles using label grid
	m_particles = new std::vector<SimUtil::Particle2D>();
	seedParticles(PARTICLES_PER_CELL, m_particles);

}

/*
Builds initial grid of dimensions (width, height) that contains the initial
geometry for the system to simulate. It reads the initial geometry from
the specified input file parameter.
Args:
width/height - grid dimensions
geomFile - the file containing the geometry
grid - the 2D array to put the initial grid in
*/
void FluidSolver2D::readInGeom(int width, int height, std::string geomFileName, Mat2Di grid) {
	// open the geometry file
	std::ifstream geomFile(geomFileName);
	if (geomFile.is_open()) {
		std::string lineStr;
		// parse file based on given dimensions, will error if file does not match these
		// fills grid so that [0][0] is at bottom left corner of simulation
		for (int i = height - 1; i >= 0; i--) {
			std::getline(geomFile, lineStr);
			for (int j = 0; j < width; j++) {
				switch (lineStr[j]) {
				case 'f':
					grid[i][j] = SimUtil::FLUID;
					break;
				case 's':
					grid[i][j] = SimUtil::SOLID;
					break;
				case 'a':
					grid[i][j] = SimUtil::AIR;
					break;
				}
			}
		}
		geomFile.close();
	}
}

/*
Seeds the initial simulation particles. Particles are created for each fluid-labeled
cell in a random-jittered subgrid pattern.
Args:
particlesPerCell - number of particles to seed in each fluid cell.
particleList - list to place the new particles in
*/
void FluidSolver2D::seedParticles(int particlesPerCell, std::vector<SimUtil::Particle2D> *particleList) {
	// set random seed
	srand(time(NULL));
	// go through all cells marked fluid
	for (int i = 0; i < m_gridHeight; i++) {
		for (int j = 0; j < m_gridWidth; j++) {
			if (m_label[i][j] == SimUtil::FLUID) {
				// seed randomly in 2x2 subgrid of the cell
				SimUtil::Vec2 cellCenter = getCellLocation(i, j);
				SimUtil::Vec2 subCenters[] = {
					SimUtil::Vec2(cellCenter.x - 0.25f*m_dx, cellCenter.y + 0.25f*m_dx), // top left
					SimUtil::Vec2(cellCenter.x + 0.25f*m_dx, cellCenter.y + 0.25f*m_dx), // top right
					SimUtil::Vec2(cellCenter.x + 0.25f*m_dx, cellCenter.y - 0.25f*m_dx), // bottom right
					SimUtil::Vec2(cellCenter.x - 0.25f*m_dx, cellCenter.y - 0.25f*m_dx) // bottom left
				};
				// cycle through subgrid to place all particles
				for (int k = 0; k < particlesPerCell; k++) {
					// randomly jitter from subgrid center
					// give a random factor from [-0.25, 0.25] multiplied by dx
					float jitterX = ((float)((rand() % 51) - 25) / 100.0f) * m_dx;
					float jitterY = ((float)((rand() % 51) - 25) / 100.0f) * m_dx;
					SimUtil::Vec2 pos(subCenters[i % 4].x + jitterX, subCenters[i % 4].y + jitterY);
					SimUtil::Vec2 vel(0.0f, 0.0f);
					particleList->push_back(SimUtil::Particle2D(pos, vel));
				}
			}
		}
	}
}

/*
Finds the physical center location of the cell with index [i][j] (ith row, jth col)
based on dx
Args:
i - row index of cell
j - col index of cell
Returns:
Vec2 (x, y) containing physical location from bottom left corner of grid.
*/
SimUtil::Vec2 FluidSolver2D::getCellLocation(int i, int j) {
	return SimUtil::Vec2(j*m_dx + 0.5f*m_dx, i*m_dx + 0.5f*m_dx);
}



