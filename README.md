# 2D Fluid Simulation
This is my first implementation of a full fluid solver. Multiple video demos of the simulation running is available [here](http://cse.unl.edu/~drempe/projects.html). I have initially done it in 2D here, but plan to extend it to 3D soon. The system is split into two components: the solver (`FluidSolver2d.cpp`) and the renderer (`FluidRenderer2d.cpp`). The solver is written from scratch based on <i>Fluid Simulation for Computer Graphics (Second Edition)</i> by Robert Bridson. Free SIGGRAPH course notes that this book is based on are [available on the author's website](https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf). The renderer mostly uses OpenGL, GLUT, and GLM in a few places. Right now it is extremely simplistic, pretty much just for debugging the solver, as my final goal is a full 3D implementation. 

![sim screenshot](https://github.com/davrempe/2d-fluid-sim/blob/master/images/screenshot_2d_fluid_sim.png "Example of simulation running")

# To Run
* The easiest way is to simply open the Visual Studio solution and run. This uses the nupengl package installed through VS.  
    * I've done all development and testing on Windows 10 with Visual Studio Community 2015, I haven't tested any other set-ups.
    * If you run into issues getting GLM to link properly, see [this Stack Overflow answer](http://stackoverflow.com/a/17912620). 
* If you don't have Visual Studio, make sure to have GLUT, GLEW, and GLM properly set up.

The simulation can be parameterized by setting variables in both `FluidSim2dMain.cpp` and `FluidSolver2d.h`. In `FluidSim2dMain.cpp` this includes MAC grid width and height (100 x 50 by default), grid cell width (0.005 m), time step size (0.01 s), input geometry, number of frames to simulate, and frame rate. In `FluidSolver2d.h` you can set the parameters for main algorithms used in the simulation like number of particles to seed per fluid cell, PIC/FLIP blending weight, CFL condition for advection, gravity, fluid density, PCG tolerance, and PCG max iterations. 

# Details
The simulation is meant to be water, as the density in the solver is set to 1000 kg/m^3 (though this can easily be changed). Therefore, the fluid is treated as inviscid and a PIC/FLIP method with an underlying MAC grid is used (again, based on the above book). At the start of the simulation, 4 particles are seeded randomly in each fluid grid cell. So in each simulation step the (high-level) process is as follows:
* transfer particle velocities to the grid
* extrapolate grid veclocity values out at least 1 cell
* save velocity grids for FLIP update
* apply body forces (gravity) on the grid
* solve for pressure on the grid
* apply pressure force on the grid
* transfer grid values to particles using a PIC/FLIP blend
* advect particles through grid velocity field

A few notes on the above steps. To transfer particle velocities to the grid, at each grid point a weighted average of nearby particles is used. "Nearby" is determined by the bilinear hat function. Body forces are applied to the grid using simple forward Euler. The pressure system is solved using the preconditioned conjugate gradient (PCG) algorithm with a Modified Incomplete Cholesky level zero, MIC(0), preconditioner. The grid to particle transfer is a mix of PIC and FLIP determined by the given weight parameter. A higher weighted PIC update tends to cause excessive diffusion, but can be used to fake more viscous fluids. Finally the particles are advected using the Runge-Kutta 3 method with adaptive substeps based on the given CFL parameter. For full details see the above mentioned course notes, and I also found [this paper](https://www.cs.ubc.ca/~rbridson/docs/zhu-siggraph05-sandfluid.pdf) very helpful. 
