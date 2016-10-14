#include <fstream>
#include <iostream>

#include "FluidSolver2d.h"
#include "SimUtil.h"


/*
Builds initial grid of dimensions (width, height) that contains the initial
geometry for the system to simulate. It reads the initial geometry from
the specified input file parameter.
Args:
gridDim - grid dimensions (width, height)
grid - the 2D array to put the initial grid in
*/
//void buildInitialGrid(Vec2 gridDim, CellType **grid) {
//	int width = (int)gridDim.x;
//	int height = (int)gridDim.y;
//
//	// open the geometry file
//	std::ifstream geomFile(INITIAL_GEOMETRY_FILE_IN);
//	if (geomFile.is_open()) {
//		std::string lineStr;
//		// parse file based on given dimensions, will error if file does not match these
//		for (int i = 0; i < height; i++) {
//			std::getline(geomFile, lineStr);
//			for (int j = 0; j < width; j++) {
//				switch (lineStr[j]) {
//				case 'f':
//					grid[i][j] = CellType::SOLID;
//					break;
//				case 's':
//
//					break;
//				case 'a':
//
//					break;
//				}
//			}
//		}
//		geomFile.close();
//	}
//}


FluidSolver2D::FluidSolver2D(int width, int height, float dx, float dt)
{
}

FluidSolver2D::~FluidSolver2D()
{
}

void FluidSolver2D::init(std::string initialGeometryFile)
{
}



