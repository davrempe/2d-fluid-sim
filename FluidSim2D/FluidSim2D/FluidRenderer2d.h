#ifndef FLUID_RENDERER_2D_H_
#define FLUID_RENDERER_2D_H_

#define BUFFER_OFFSET(offset) ((GLvoid*) offset)

#include <gl/glew.h>
#include <glm/glm.hpp>
#include <string>
#include <vector>

#include "SimUtil.h"

class FluidRenderer2D {

	struct VertexData {
		GLfloat vPos[2];
		GLfloat vColor[3];

		VertexData() {
			vPos[0] = 0.0;
			vPos[1] = 0.0;
			vColor[0] = 0.0;
			vColor[1] = 0.0;
			vColor[2] = 0.0;
		}

		VertexData(GLfloat vertPos[2], GLfloat vertColor[3]) {
			vPos[0] = vertPos[0];
			vPos[1] = vertPos[1];
			vColor[0] = vertColor[0];
			vColor[1] = vertColor[1];
			vColor[2] = vertColor[2];
		}
	};

private:
	//----------------------------------------------------------------------
	// Rendering Parameters
	//----------------------------------------------------------------------

	// file name of geometry information
	std::string m_geomFile;
	// file name of particle data from simulation
	std::string m_particleFile;
	// width of simulation grid
	int m_width;
	// height of simulation grid
	int m_height;
	// grid cell dimensions
	float m_dx;

	//----------------------------------------------------------------------
	// Rendering-related members
	//----------------------------------------------------------------------

	GLfloat PARTICLE_COLOR[3] = { 0.0f, 0.0f, 1.0f };
	GLfloat SOLID_COLOR[3] = { 1.0f, 1.0f, 1.0f };
	GLfloat BACKGROUND_COLOR[3] = { 0.2f, 0.2f, 0.2f };

	const int VERTICES_PER_QUAD = 4;
	const int INDICES_PER_QUAD = 6;
	const int BYTES_PER_FLOAT = sizeof(GLfloat);

	// frame currently being rendered
	int m_currentFrame;
	// the seconds between each frame for the given frame rate
	float m_frameTime;

	// particle position data for every frame
	std::vector<std::vector<SimUtil::Vec2>> m_particlePosData;
	// grid of labels representing solid geometry
	SimUtil::Mat2Di m_geomGrid;
	// transform matrix
	glm::mat4 m_transMat;

	// array of VertexData for all particles in current frame
	VertexData *m_particleVertData;
	// number of particles in current particle VertexData array
	int m_numParticlesInFrame;
	// array of VertexData for vertices that make up solids
	VertexData *m_solidVertData;
	// array of indices for how to connect solid vertices
	GLushort *m_solidIndData;
	// number of solid quads in current vertex/index data arrays
	int m_numberSolidCells;
	// offset of color attribute in VertexData struct
	GLintptr m_vColorOffset;

	// buffer for particle vertices
	GLuint m_pvBuffer;
	// buffer for solid vertices
	GLuint m_svBuffer;
	// buffer for solid indices
	GLuint m_siBuffer;

	// shader position attribute
	GLuint m_vPos;
	// shader color attribute
	GLuint m_vColor;
	// uniform shader projection matrix attribute
	GLuint m_MVP;
	// shader program
	GLuint m_shaderProgram;

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	void readInParticleData();
	void strSplit(const std::string&, char, std::vector<std::string>&);
	void updateParticleVertexData(int);
	void initSolidVertexData();
	void updateBufferData(int);
	void initGL();
	// time function to control frame rate
	// TODO
	
public:
	/*
	Creates a new 2D fluid renderer.
	Args:
	- geomFile - name of the file with geometry data in it. Used to draw solid objects in scene.
	- particleFile - name of the file with particle data in it.
	- frameRate - the frame rate to render at
	- gridWidth - the width of the simulation grid used
	- gridHeight - the height of the simulation grid used
	- cellWidth - width of a cell in the simulation grid used
	*/
	FluidRenderer2D(std::string, std::string, float, int, int, float);
	~FluidRenderer2D();

	/*
	Initializes renderer by reading in geometry/particle data and preparing to render.
	Args:
	argc - unmodified argc from main
	argv - unmodified argv from main
	*/
	void init(int, char**);

	/*
	Starts the renderer.
	*/
	void render();

	// display function, should not be called from outside the renderer
	void display();
	// time function, should not be called from outside the renderer
	void timer(int);
};

#endif //FLUID_RENDERER_2D_H_