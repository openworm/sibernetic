Smoothed-Particle-Hydrodynamics
===============================

This is a C++ implementation of the [Smoothed Particle Hydrodynamics](http://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) algorithm for the [OpenWorm](http://openworm.org) project.

There is also a separate effort to port this code to Java, as part of the Geppetto simulation framework. 

This branch includes development on electrophysiology integration. It should currently be considered unstable.

For record configuraton into file you need to run simulation with key -l_to - it create 3 new files 
at the folder ./buffers:
- connection_buffers.txt - it need to store information about conection among of elastic partciles
- membranes_buffer.txt   - it need to store information about membranes 
- position_buffer.txt    - it need to store information current position all of the non boundary particles it save information to this file every 10 steps of simulation. You shoulld remember that than more info you 
want to store than bigger output file is. 

For view result you should run simulation with key -l_from - it get positions from position_buffer.txt file and 
draw evolution of system in time
