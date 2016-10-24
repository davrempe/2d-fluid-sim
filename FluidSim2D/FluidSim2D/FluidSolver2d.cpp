#include <fstream>
#include <iostream>

#include "FluidSolver2d.h"
#include "SimUtil.h"


FluidSolver2D::FluidSolver2D(int width, int height, float dx, float dt){
	this->gridWidth = width;
	this->gridHeight = height;
	this->dx = dx;
	this->dt = dt;
}

FluidSolver2D::~FluidSolver2D(){
	// clean up grid
	SimUtil::deleteMat2D<int>(gridHeight, gridWidth, label);
	SimUtil::deleteMat2D<float>(gridHeight, gridWidth, p);
	SimUtil::deleteMat2D<float>(gridHeight, gridWidth + 1, u);
	SimUtil::deleteMat2D<float>(gridHeight + 1, gridWidth, v);
}

void FluidSolver2D::init(std::string initialGeometryFile){
	// set up the grid for simulation by initializing arrays
	label = SimUtil::initMat2D<int>(gridHeight, gridWidth);
	p = SimUtil::initMat2D<float>(gridHeight, gridWidth);
	u = SimUtil::initMat2D<float>(gridHeight, gridWidth + 1);
	v = SimUtil::initMat2D<float>(gridHeight + 1, gridWidth);

	// read in initial geometry to populate label grid
	readInGeom(gridWidth, gridHeight, initialGeometryFile, label);
	//SimUtil::printMat2D<int>(gridHeight, gridWidth, label);


}

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



