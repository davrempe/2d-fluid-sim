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
// the number of frames to simulate
const int NUM_SIM_FRAMES = 125;

//----------------------------------------------------------------------
// I/O Parameters
//----------------------------------------------------------------------

// input file for initial system state - grid marked solid, fluid, or air
const std::string INITIAL_GEOMETRY_FILE_IN = "initial_geometry.txt";
// output file for particle data
const std::string PARTICLE_DATA_FILE_OUT = "particle_data.txt";

//----------------------------------------------------------------------
// Global Variables
//----------------------------------------------------------------------

int main() {
	if (RUN_SIM) {
		FluidSolver2D *solver = new FluidSolver2D(GRID_WIDTH, GRID_HEIGHT, GRID_CELL_WIDTH, TIME_STEP);
		

		// cleanup
		delete solver;
	}

	if (RUN_RENDERING) {

	}


	return 0;
}