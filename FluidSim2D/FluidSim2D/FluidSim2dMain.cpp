//----------------------------------------------------------------------
// Title: 2D Fluid Simulation
// Author: Davis Rempe
//
// Implementation based on "Fluid Simulation for Computer Graphics"
// by Robert Bridson.
//----------------------------------------------------------------------

#include <string>
#include <fstream>
#include <iostream>

#include "FluidSolver2d.h"
#include "FluidRenderer2d.h"
#include "SimUtil.h"

//----------------------------------------------------------------------
// Execution Options
//----------------------------------------------------------------------

// whether to run the simulation
const bool RUN_SIM = true;
// whether to run rendering
const bool RUN_RENDERING = true;

//----------------------------------------------------------------------
// Simulation Parameters
//----------------------------------------------------------------------

// resolution of the grid to use (width, height)
const int GRID_WIDTH = 100;
const int GRID_HEIGHT = 50;
// grid cell width (in meters)
const float GRID_CELL_WIDTH = 0.005f;
// simulation time step (in seconds)
const float TIME_STEP = 0.04f;

//----------------------------------------------------------------------
// I/O Parameters
//----------------------------------------------------------------------

// input file for initial system state - grid marked solid, fluid, or air
const std::string INITIAL_GEOMETRY_FILE_IN = "initial_geometry.txt";
// output file for particle data
const std::string PARTICLE_DATA_FILE_OUT = "particle_data.csv";
// the number of frames to simulate
const int NUM_SIM_FRAMES = 125;
// frame rate for render (fps)
const float FRAME_RATE = 25.0f;
// time step between outputted frames
const float FRAME_TIME_STEP = 1.0f / FRAME_RATE;

//----------------------------------------------------------------------
// Global Variables
//----------------------------------------------------------------------



int main() {
	if (RUN_SIM) {
		// open and clear output file
		std::ofstream *particleOut = new std::ofstream(PARTICLE_DATA_FILE_OUT, std::ofstream::trunc);

		FluidSolver2D solver(GRID_WIDTH, GRID_HEIGHT, GRID_CELL_WIDTH, TIME_STEP);
		solver.init(INITIAL_GEOMETRY_FILE_IN);
		
		// run simulation
		int framesOut = 0;
		float t = 0;
		while (framesOut < NUM_SIM_FRAMES) {
			// perform sim time step
			solver.step();

			t += TIME_STEP;
			// check if need to output data to render at current time
			if (t - ((framesOut + 1) * FRAME_TIME_STEP) <= 0.000001f) {
				solver.saveParticleData(particleOut);
				framesOut++;
			}
		}

		// cleanup
		particleOut->close();
		delete particleOut;
	}

	if (RUN_RENDERING) {
		FluidRenderer2D renderer(INITIAL_GEOMETRY_FILE_IN, PARTICLE_DATA_FILE_OUT, GRID_WIDTH, GRID_HEIGHT, GRID_CELL_WIDTH);
		renderer.init();

	}


	return 0;
}