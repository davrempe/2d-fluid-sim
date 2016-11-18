#include "FluidRenderer2d.h"
#include "FluidSolver2d.h"
#include "LoadProgram.h"

#include <gl/glut.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

// current instance of the renderer
FluidRenderer2D* g_rendererContext;

extern "C"
void displayCallback() {
	g_rendererContext->display();
}

extern "C"
void timerCallback(int id) {
	g_rendererContext->timer(id);
}

FluidRenderer2D::FluidRenderer2D(std::string geomFile, std::string particleFile, float frameRate, int gridWidth, int gridHeight, float cellWidth) : m_geomFile(geomFile), m_particleFile(particleFile), m_width(gridWidth), m_height(gridHeight), m_dx(cellWidth) {
	// initialize vertex data arrays
	m_particleVertData = new VertexData[1];

	m_currentFrame = 0;
	m_frameTime = 1.0f / frameRate;
}

FluidRenderer2D::~FluidRenderer2D() {
	delete[] m_particleVertData;
	delete[] m_solidVertData;
	delete[] m_solidIndData;
}

void FluidRenderer2D::init(int argc, char** argv) {
	// read in particle data to array
	readInParticleData();
	// read in geometry data
	m_geomGrid = SimUtil::initGrid2D<int>(m_width, m_height);
	SimUtil::readInGeom(m_width, m_height, m_geomFile, m_geomGrid);

	// set up glut and glew
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(720, 720);
	glutCreateWindow("Fluid Simulation");
	glewInit();

	// set up OpenGL
	initGL();

	g_rendererContext = this;
	glutDisplayFunc(displayCallback);
	/*glutReshapeFunc(reshape);
	glutKeyboardFunc(key);
	glutSpecialFunc(specialKey);
	glutMouseFunc(mouse);
	glutMotionFunc(drag);
	glutPassiveMotionFunc(trackMousePos);*/

}

void FluidRenderer2D::render() {
	glutMainLoop();
}

//----------------------------------------------------------------------
// Private Functions
//----------------------------------------------------------------------

/*
Fills array of particle VertexData with particle data from the given frame.
If the frame is greater than the total number of frames, it will show the
given frameNum % totalNumFrames.
Args:
- frameNum - the frame to update the array with
*/
void FluidRenderer2D::updateParticleVertexData(int frameNum) {
	int totalFrames = m_particlePosData.size();
	// clear array from previous frame
	delete[] m_particleVertData;
	// reallocate to proper size
	m_numParticlesInFrame = m_particlePosData[frameNum % totalFrames].size();
	m_particleVertData = new VertexData[m_numParticlesInFrame];

	// populate array
	for (int i = 0; i < m_numParticlesInFrame; i++) {
		SimUtil::Vec2 curParticle = m_particlePosData[frameNum % totalFrames][i];
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
	for (int i = 0; i < m_width; i++) {
		for (int j = 0; j < m_height; j++) {
			if (m_geomGrid[i][j] == SimUtil::SOLID) {
				// build quad over that cell
				SimUtil::Vec2 cellCenter = SimUtil::getGridCellPosition(i, j, m_dx);
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
Updates the buffers for the current frame being rendered.
*/
void FluidRenderer2D::updateBufferData(int frameNum) {
	// first update the particle data for this frame
	updateParticleVertexData(frameNum);

	// refill the buffer with new data
	glBindBuffer(GL_ARRAY_BUFFER, m_pvBuffer);
	glBufferSubData(GL_ARRAY_BUFFER, 0, m_numParticlesInFrame * sizeof(VertexData), m_particleVertData);

	// solid data does not change frame to frame
}

/*
Initializes OpenGL for use by initializing buffers, shaders, and attributes.
*/
void FluidRenderer2D::initGL() {
	// enable blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPointSize(4);

	// fill particle vertex array with initial data
	updateParticleVertexData(0);
	// fill solids vertex array with initial data
	initSolidVertexData();
	// create transformation matrix
	float scale = 2.0f / (m_width*m_dx);
	glm::mat4 scaleMat = glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, 1.0f));
	glm::mat4 translateMat = glm::translate(glm::mat4(1.0f), glm::vec3((-m_width*m_dx) / 2.0f, (-m_height*m_dx) / 2.0f, 0.0f));
	m_transMat = scaleMat* translateMat;

	// calculate offset for color in VertexData object
	m_vColorOffset = BYTES_PER_FLOAT * 2;

	// create buffers
	GLuint buffs[3];
	glGenBuffers(3, buffs);
	m_pvBuffer = buffs[0];
	m_svBuffer = buffs[1];
	m_siBuffer = buffs[2];

	//init buffers
	// particle vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, m_pvBuffer);
	glBufferData(GL_ARRAY_BUFFER, m_numParticlesInFrame * sizeof(VertexData), m_particleVertData, GL_DYNAMIC_DRAW);
	// solid vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, m_svBuffer);
	glBufferData(GL_ARRAY_BUFFER, m_numberSolidCells * VERTICES_PER_QUAD * sizeof(VertexData), m_solidVertData, GL_STATIC_DRAW);
	// solid index buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_siBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_numberSolidCells * INDICES_PER_QUAD * sizeof(GLushort), m_solidIndData, GL_STATIC_DRAW);

	// create shaders
	const char* vShadeFile = "vert.glsl";
	const char* fShadeFile = "frag.glsl";

	m_shaderProgram = LoadProgram(vShadeFile, fShadeFile);
	m_vPos = glGetAttribLocation(m_shaderProgram, "vPos");
	m_vColor = glGetAttribLocation(m_shaderProgram, "vColor");
	m_MVP = glGetUniformLocation(m_shaderProgram, "MVP");	

	glClearColor(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2], 1.0f);
}

// called at certain time increments to redraw and calculate fps
void FluidRenderer2D::timer(int id) {
	glutPostRedisplay();
}

/*
Renders current particle and solid data.
*/
void FluidRenderer2D::display() {
	
	// update buffer data
	updateBufferData(m_currentFrame);

	glClear(GL_COLOR_BUFFER_BIT);

	glUseProgram(m_shaderProgram);
	// draw solids
	glBindBuffer(GL_ARRAY_BUFFER, m_svBuffer);
	glVertexAttribPointer(m_vPos, 2, GL_FLOAT, GL_FALSE, sizeof(VertexData), BUFFER_OFFSET(0));
	glVertexAttribPointer(m_vColor, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), BUFFER_OFFSET(m_vColorOffset));
	glEnableVertexAttribArray(m_vPos);
	glEnableVertexAttribArray(m_vColor);
	glUniformMatrix4fv(m_MVP, 1, GL_FALSE, glm::value_ptr(m_transMat));

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_siBuffer);
	glDrawElements(GL_TRIANGLES, m_numberSolidCells * INDICES_PER_QUAD, GL_UNSIGNED_SHORT, BUFFER_OFFSET(0));

	// draw particles
	glBindBuffer(GL_ARRAY_BUFFER, m_pvBuffer);
	glVertexAttribPointer(m_vPos, 2, GL_FLOAT, GL_FALSE, sizeof(VertexData), BUFFER_OFFSET(0));
	glVertexAttribPointer(m_vColor, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), BUFFER_OFFSET(m_vColorOffset));
	glEnableVertexAttribArray(m_vPos);
	glEnableVertexAttribArray(m_vColor);
	glUniformMatrix4fv(m_MVP, 1, GL_FALSE, glm::value_ptr(m_transMat));

	glDrawArrays(GL_POINTS, 0, m_numParticlesInFrame);

	// free buffers
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// actually draw
	glutSwapBuffers();

	m_currentFrame++;

	// call timer function
	// set redisplay at inputted fps
	glutTimerFunc(m_frameTime*1000, timerCallback, 0);
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

			if (frameVec.size() > 0) {
				m_particlePosData.push_back(frameVec);
			}
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