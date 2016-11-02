#include "FluidRenderer2d.h"
#include "FluidSolver2d.h"

#include <gl/glut.h> 
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

FluidRenderer2D::FluidRenderer2D(std::string geomFile, std::string particleFile, int gridWidth, int gridHeight, float cellWidth) : m_geomFile(geomFile), m_particleFile(particleFile), m_width(gridWidth), m_height(gridHeight), m_dx(cellWidth) {
	// initialize vertex data arrays
	m_particleVertData = new VertexData[1];
}

FluidRenderer2D::~FluidRenderer2D() {
	delete[] m_particleVertData;
	delete[] m_solidVertData;
	delete[] m_solidIndData;
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
	// clear array from previous frame
	delete[] m_particleVertData;
	// reallocate to proper size
	m_numParticlesInFrame = m_particlePosData[frameNum].size();
	m_particleVertData = new VertexData[m_numParticlesInFrame];

	// populate array
	for (int i = 0; i < m_numParticlesInFrame; i++) {
		SimUtil::Vec2 curParticle = m_particlePosData[frameNum][i];
		GLfloat particlePos[2] = { curParticle.x, curParticle.y };
		m_particleVertData[i] = VertexData(particlePos, PARTICLE_COLOR);
	}
}

/*
Fills array of solid VertexData and index data with the solid objects.
*/
void FluidRenderer2D::initSolidVertexData() {
	// max size is if every cell in the grid solid
	m_solidVertData = new VertexData[m_width*m_height * VERTICES_PER_QUAD];
	m_solidIndData = new GLushort[m_width*m_height * INDICES_PER_QUAD];

	// traverse geometry grid to find solid cells
	m_numberSolidCells = 0;
	for (int i = 0; i < m_height; i++) {
		for (int j = 0; j < m_width; j++) {
			if (m_geomGrid[i][j] == SimUtil::SOLID) {
				// build quad over that cell
				SimUtil::Vec2 cellCenter = SimUtil::getCellLocation(i, j, m_dx);
				GLfloat nwPos[2] = { cellCenter.x - 0.5f*m_dx, cellCenter.y + 0.5f*m_dx };
				GLfloat nePos[2] = { cellCenter.x + 0.5f*m_dx, cellCenter.y + 0.5f*m_dx };
				GLfloat sePos[2] = { cellCenter.x + 0.5f*m_dx, cellCenter.y - 0.5f*m_dx };
				GLfloat swPos[2] = { cellCenter.x - 0.5f*m_dx, cellCenter.y - 0.5f*m_dx };
				// add to vert array
				GLushort nwInd = m_numberSolidCells*VERTICES_PER_QUAD + 0;
				GLushort neInd = m_numberSolidCells*VERTICES_PER_QUAD + 1;
				GLushort seInd = m_numberSolidCells*VERTICES_PER_QUAD + 2;
				GLushort swInd = m_numberSolidCells*VERTICES_PER_QUAD + 3;
				m_solidVertData[nwInd] = VertexData(nwPos, SOLID_COLOR);
				m_solidVertData[neInd] = VertexData(nePos, SOLID_COLOR);
				m_solidVertData[seInd] = VertexData(sePos, SOLID_COLOR);
				m_solidVertData[swInd] = VertexData(swPos, SOLID_COLOR);
				// add to index array
				m_solidIndData[m_numberSolidCells*INDICES_PER_QUAD + 0] = nwInd;
				m_solidIndData[m_numberSolidCells*INDICES_PER_QUAD + 1] = swInd;
				m_solidIndData[m_numberSolidCells*INDICES_PER_QUAD + 2] = seInd;
				m_solidIndData[m_numberSolidCells*INDICES_PER_QUAD + 3] = nwInd;
				m_solidIndData[m_numberSolidCells*INDICES_PER_QUAD + 4] = seInd;
				m_solidIndData[m_numberSolidCells*INDICES_PER_QUAD + 5] = neInd;

				m_numberSolidCells++;
			}
		}
	}
}

/*
Initializes OpenGL for use by initializing buffers, shaders, and attributes.
*/
void FluidRenderer2D::initGL() {
	// enable blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// fill particle vertex array with initial data
	updateParticleVertexData(0);
	// fill solids vertex array with initial data
	initSolidVertexData();
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