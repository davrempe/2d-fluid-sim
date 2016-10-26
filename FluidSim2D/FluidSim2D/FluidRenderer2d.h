#ifndef FLUID_RENDERER_2D_H_
#define FLUID_RENDERER_2D_H_

#include <gl/glew.h>

class FluidRenderer2D {

	struct VertexData {
		GLfloat vPos[3];
		GLfloat vColor[3];

		VertexData() {
			vPos[0] = 0.0;		vPos[1] = 0.0;		vPos[2] = 0.0;
			vColor[0] = 0.0;		vColor[1] = 0.0;		vColor[2] = 0.0;
		}

		VertexData(GLfloat vertPos[3], GLfloat vertColor[3]) {
			vPos[0] = vertPos[0];		vPos[1] = vertPos[1];		vPos[2] = vertPos[2];
			vColor[0] = vertColor[0];		vColor[1] = vertColor[1];		vColor[2] = vertColor[2];
		}
	};

private:
	// array holding particle position data for each frame

	// array of VertexData for all particles in current frame
	// array of VertexData for vertices that make up solids
	// array of indices for how to connect solid vertices

	// buffer for particle vertices
	// buffer for solid vertices
	// buffer for solid indices

	// shader position attribute
	// shader color attribute

	// offset for color within VertexData struct

	// function to read in particle position data from file, populate array
	// function to update particle VertexData array for GL each frame (should just be looking at next frame's position data and translating to GL)
	// function to actually update buffer with new particle data
	// init openGL function
	// time function to control frame rate
	// display function

	// render function (basically like main)

public:
	FluidRenderer2D();
	~FluidRenderer2D();


};

#endif //FLUID_RENDERER_2D_H_