#include "FluidRenderer2d.h"
#include "FluidSolver2d.h"

#include <gl/glut.h> 
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

FluidRenderer2D::FluidRenderer2D(std::string geomFile, std::string particleFile, int gridWidth, int gridHeight, float cellWidth) : m_geomFile(geomFile), m_particleFile(particleFile), m_width(gridWidth), m_height(gridHeight), m_dx(cellWidth) {
	
}

FluidRenderer2D::~FluidRenderer2D() {

}

void FluidRenderer2D::init() {
	// read in particle data to array
	readInParticleData();
	// read in geometry data
	m_geomGrid = SimUtil::initMat2D<int>(m_height, m_width);
	SimUtil::readInGeom(m_width, m_height, m_geomFile, m_geomGrid);

	// set up OpenGL
	initGL();
}

//----------------------------------------------------------------------
// Private Functions
//----------------------------------------------------------------------

/*
Fills array of particle VertexData with particle data from the given frame.
Args:
- frameNum - the frame to update the array with
*/
void FluidRenderer2D::updateParticleVertexData(int frameNum) {

}

/*
Fills array of solid VertexData and index data with the solid objects.
*/
void FluidRenderer2D::updateSolidVertexData() {

}

/*
Initializes OpenGL for use by initializing buffers, shaders, and attributes.
*/
void FluidRenderer2D::initGL() {
	// fill particle vertex array with initial data

	// fill solids vertex array with initial data

	// create projection matrix
}

/*
Reads in particle position data from the particle file and uses it to initialize
the particle position data array.
*/
void FluidRenderer2D::readInParticleData() {
	std::cout << "Reading in particle data for..." << std::endl;
	// open the particle data file
	std::ifstream particleFile(m_particleFile);
	if (particleFile.is_open()) {
		int curFrame = 0;
		while (particleFile.good()) {
			std::cout << "Frame " << curFrame << std::endl;
			std::vector<SimUtil::Vec2> frameVec;
			// loop through all frames and fill particle data array
			std::string lineStr;

			std::getline(particleFile, lineStr);
			std::vector<std::string> tokens;
			strSplit(lineStr, ',', tokens);
			for (int i = 0; i < tokens.size(); i++) {
				// parse particle position (e.g. "0.0345 0.1235")
				std::vector<std::string> strPos;
				strSplit(tokens[i], ' ', strPos);
				SimUtil::Vec2 pos(atof(strPos[0].c_str()), atof(strPos[1].c_str()));
				frameVec.push_back(pos);
			}

			m_particlePosData.push_back(frameVec);
			curFrame++;
		}
		particleFile.close();
	}
}

/*
Splits a string based on delimeter.
Grabbed from http://stackoverflow.com/a/236803.
*/
void FluidRenderer2D::strSplit(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
}