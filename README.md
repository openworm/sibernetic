![Sibernetic](http://i.imgur.com/Hbsw6Zs.png)

Sibernetic is physical simulator of biomechanical matter (membranes, elastic matter, contractile matter) and environments (liquids, solids and elastic matter with variable physical properties) developed for simulations of C. elegans physical body dynamics within the [OpenWorm project](http://www.openworm.org) by Andrey Palyanov, Sergey Khayrulin and Mike Vella (development of a Python module for external muscle activating signals generation and input) as part of the [OpenWorm team](http://www.openworm.org/people.html). At its core, Sibernetic is built as an extension to Predictive-Corrective Incompressible Smoothed Particle Hydrodynamics (PCISPH). It is primarily written in  C++ and OpenCL, which makes possible to run simulations on CPUs or GPUs, and has 3D visualization support built on top of OpenGL.

There is a separate effort lead by [Giovanni Idili](https://github.com/gidili) and [Sergey Khayrulin](https://github.com/skhayrulin) to port this code to Java, as part of the [Geppetto simulation framework](http://www.geppetto.org/).

Compiling / running (Linux/mac)
------------------------------

[![Join the chat at https://gitter.im/openworm/sibernetic](https://badges.gitter.im/openworm/sibernetic.svg)](https://gitter.im/openworm/sibernetic?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

**Linux**
Install OpenCL on Ubuntu. We suggest you initially go with [AMD OpenCL drivers](http://developer.amd.com/tools-and-sdks/heterogeneous-computing/amd-accelerated-parallel-processing-app-sdk/downloads/) as we have found these to be the most stable and complete. You can also try [Intel's drivers](http://develnoter.blogspot.co.uk/2012/05/installing-opencl-in-ubuntu-1204.html). This step often causes problems, contact the [openworm-discuss](https://groups.google.com/forum/#!forum/openworm-discuss) mailing list if you encounter issues. The AMD drivers include samples in `/opt/AMDAPP/samples/opencl/bin` which you can use to verify your OpenCL support is working.

You'll also need a variety of libraries. In ubuntu, install the dependencies with:

```bash
sudo apt-get install g++ python-dev freeglut3-dev nvidia-opencl-dev libglu1-mesa-dev libglew-dev python-numpy
```

Next, from the `sibernetic/` directory run:

```bash
make clean
make all
```

Also you  may need to set some enviromat variables like path to OpenCL lib or header for to do this
you can fix your `LD_LIBRARY_PATH` as this:

```bash
export LD_LIBRARY_PATH=/path/to/opencl_lib/folder/:$LD_LIBRARY_PATH
e.g.
export LD_LIBRARY_PATH=/usr/local/cuda/lib64/:$LD_LIBRARY_PATH
```

You can find OpenCL lib in CUDA folder if you're using NVIDIA (`/usr/local/cuda/lib64/`) or run this command.

```bash
ldconfig -p | grep opencl
```

Also you may need to give compiler path to OpenCL header files usually you can find them in `/usr/include/CL` if they there than you don't need do anything. In othe case you can edit makefile directly and add directory to OpenCL headers by adding options `-I/path/to/opencl_includes/` or you can copy folder with header into `/usr/include/` but you should have root permission for doing that.

**Mac**: stay in the top-level folder. You need before run export several environment variables:

```bash
export PYTHONHEADERDIR=/usr/local/Cellar/python/<version_of_installed_pythonFramework>/Python.framework/Headers/
export PYTHONLIBDIR=/usr/local/lib/python2...
export PYTHONFRAMEWORKDIR=/usr/local/Frameworks/
```

Then

```bash
make clean -f makefile.OSX
make all -f makefile.OSX
```

You should see an output which looks something like this:

```bash
Building file: ../src/PyramidalSimulation.cpp
Invoking: GCC C++ Compiler

....
more stuff...
....

Building target: Sibernetic
Invoking: GCC C++ Linker
g++ -L/usr/lib -L/usr/lib/python2.7 -o "Sibernetic"  ./src/PyramidalSimulation.o ./src/main.o ./src/owHelper.o ./src/owOpenCLSolver.o ./src/owPhysicsFluidSimulator.o ./src/owWorldSimulation.o   -lOpenCL -lpython2.7 -lrt -lglut -lGL -lGLU
Finished building target:Sibernetic
```

Then navigate to the top-level folder in the hierarchy (e.g `Sibernetic`) and set your `PYTHONPATH`:

```bash
export PYTHONPATH=$PYTHONPATH:.
```

Finally, to run, run the command:

**Linux**:

```bash
./Release/Sibernetic
```

**Mac**:

```bash
./Release/Sibernetic
```

You may need to make `./Release/Sibernetic` executable like so:

```bash
chmod +x ./Release/Sibernetic
```

If you do not run from the top-level folder you will see an error which looks something like this:

```bash
Compilation failed:
"/tmp/OCLQ1BaOw.cl", line 8: catastrophic error: cannot open source file
"src//owOpenCLConstant.h"
#include "src//owOpenCLConstant.h"
```

What's inside
-------------

Physical Algorithms:

- PCI SPH - simulation incompressible liquid [1]
- Simulation elastic matter
- Simulation liquid-impermeable membranes
- Boundary handling [2]
- Surface tension [3]

There are two demo scenes generated for Sibernetic. The first one contains an elastic cube covered with liquid-impermeable membranes and liquid inside. The second one contains two elastic membranes attached to a boundary (one of them has liquid-impermeable membranes covering them and the other one doesn't).

To switch between demos you need to press the 1 or 2 keys respectively. To pause simulation you may press space bar.

References

1. B. Solenthaler, Predictive-Corrective Incompressible SPH. ACM Transactions on Graphics (Proceedings of SIGGRAPH), 28(3), 2009.
2. M. Ihmsen, N. Akinci, M. Gissler, M. Teschner, Boundary Handling and Adaptive Time-stepping for PCISPH Proc. VRIPHYS, Copenhagen, Denmark, pp. 79-88, Nov 11-12, 2010.
3. M. Becker, M. Teschner. Weakly compressible SPH for free surface flows // Proceedings of the 2007 ACM SIGGRAPH/Eurographics symposium on Computer animation, pages 209-217.

Main command options
--------------
To start Sibernetic with argument print in command prompt next `./Release/Sibernetic -whatever
Available options:
```
 -no_g                 Run without graphics
 -l_to                 Save simulation results to disk.
 -export_vtk           Save simulation results to VTK files.
     logstep=<value>   Log every <value> steps
 -l_from               Load simulation results from disk.
     lpath=<value>     Indicates path to the directory (not the file) where result of simulation will be stored.  
                       This option work only for -l_to and -l_from options
 -test                 Run some physical tests.
 -f <filename>         Load configuration from file <filename>.
 device=<device_type>  Trying to init OpenCL on device <type> it could be cpu or gpu
                       default-ALL (it will try to init most powerful available device).
 timestep=<value>      Start simulation with time step = <value> in seconds.
 timelimit=<value>     Run simulation until <value> will be reached in seconds.
 leapfrog              Use for integration LeapFrog method
 oclsourcepath=<value> You can indicate path to you'r OpenCL program just using this option
 -nrn <value>          Indicates that you plan run simulation with NEURON simulation = <value> value should be a file which
                       can be run by NEURON simulator and also you should have installed neuron and sibernetic_neuron bridge.
 -help                 Print this information on screen.
```

LeapFrog integration
--------------
[Leapfrog](https://en.wikipedia.org/wiki/Leapfrog_integration) is second order method insted of [Semi-implicid Euler](https://en.wikipedia.org/wiki/Semi-implicit_Euler_method) which we are using as default method for integration. For run simulation with Leapfro integration medhod print run command
```
./Release/Sibernetic leapfrog
```

Run simulation from configuration file
--------------
All configuration is stored in the [configuration folder](configuration). There are two demo configurations [demo1](configuration/demo1)
and [demo2](configuration/demo2) (demo1 is the default configuration). You can switch between two demo configurations
directly inside the working Sibernetic - just push button '1' or '2' respectively. To run your configuration put your configuration file
into the configuration folder and run Sibernetic using:
```
./Release/Sibernetic -f <configuration_file_name>.
```
To run the worm body simulation you need run Sibernetic with key:
```
./Release/Sibernetic -f worm
```
It loads the [worm body configuration](configuration/worm) and initialises and runs the [Python module](src/main_sim.py) which is
responsible for muscle signal updating.
For run simulation with crawling worm on carpet like surface or swimming in deep water you need run sibenretic with next command arguments:
```
./Release/Sibernetic -f worm_crawling oclsourcepath=src/sphFluid_crawling.cl
./Release/Sibernetic -f worm_deep_water oclsourcepath=src/sphFluid_crawling.cl
```

Control in graphical mode
---------------
If you run Sibernetic with graphics you can work with scene rotation and scaling using the mouse. There are also several control button options available:
```
'Space' - pause the simulation
's'     - save current configuration into file
          ./configuration/snapshot/configuration_name_current_time_and_date you can run this
than (./Release/Sibernetic -f ./configuration/snapshot/configuration_default).
'q' or 'Esc'     - quit the sibernetic
'1'     - run demo1 configuration
'2'     - run demo2 configuration
'r'     - reset current configuration and start from begining
```

Configuration file format
---------------
The configuration file consists of:
First block is an optional if you didn't indicate this block then sibernetic will init consts by default value which you can find in [owPhysicsConstant.h ](./inc/owPhysicsConstant.h).
```
[physical parameters]
mass: 5.4e-14
timeStep: 5.0e-06
simulationScale: 2.46e-06
viscosity: 5.0e-05
surfTensCoeff: 1.21948e+27
elasticityCoefficient: 5.55556e+08
```
Next 6 lines is a spatial description of boundary box
```
[simulation box]
xmin
xmax
ymin
ymax
zmin
zmax
[position] - contains information about position of all particles e.g.
1 0 0 1
1 0 1 1
...
[velocity] - contains information about velocities of all particles e.g.
0 0 0 1
0 0 0 1
...
[connection] - contains information about elastic connection of all elastic particles e.g.
1	1.58649939377	1.1	0.0 
7	1.58649939377	1.1	0.0
...
[membranes] - contains information about membranes e.g.
0	1	7
7	8	1
...
[particleMemIndex] - contains information about in which membranes elastic particle is includes e.g.
0
144
288
-1
...
```
Position and velocity are represented as 4D vectors which contains information about x, y, z of particle in space and information about particle's type (it could be liquid - 1, elastic - 2 or boundary - 3). Each elastic particle has 32 places allocated in connections buffer. Each connection is represented like a 4D vector 
ID of particle to connected to 
stedy-state lenght of connection 
id of muscle if this connection is a muscle fiber
and unused data - need for vectorization.
Connections buffer stored in memory like 1D vector: length of each is equal to `NUM_OF_ELASTIC_PARTICLES * 32 * 4`. So for each particular elastic particle you can find information for elastic its connections simply get sub-buffer of connection from `INTRESTING_PARTICLE_ID * 32 * 4` to `INTRESTING_PARTICLE_ID * 32 * 4 + 32 * 4`. 
Each membrane is defined by 3 elastic particles and contains 3 IDs of this particles. particleMemIndex - contains IDs of membrane in which each elastic particle is included we suppose that max numbers of membrane for one particle is 7 so particleMemIndex contains `7 * NUM_OF_ELASTIC_PARTICLES` and you can get interesting information from this buffer just get sub-buffer from indexes `INTRESTING_PARTICLE_ID * 7` to `INTRESTING_PARTICLE_ID * 7 + 7`.

Saving to disk
--------------
You can run Sibernetic on GPU. For this you should start Sibernetic with key:
```
./Release/Sibernetic device=gpu
```

You may wish to save simulations to disk rather than visualise them (**WARNING**: This is buggy)

To record configurations to file you need to run simulation with key -l_to:
```
./Release/Sibernetic -l_to
```

This create 3 new files in the folder `./buffers`:
- connection_buffers.txt - stores information about connection among the elastic particles
- membranes_buffer.txt   - stores information about membranes
- position_buffer.txt    - stores information about current position of all of the non boundary particles it save information to this file every 10 steps of simulation. You should remember that the more info you
- pressure_buffer.txt    - stores information oabout pressure for all shell particles.
want to store than bigger output file is.

For view result you should run simulation with:
```
./Release/Sibernetic -l_from
```
It get positions from position_buffer.txt file and displays the evolution of system in time


Output to the VTK files
-----------------------
Results of the simulation can be saved to the VTK files, allowing visualisation e.g. in [Paraview](http://www.paraview.org/).
To save the VTK files, run the program with the parameter `-export_vtk`,
```
./Release/Sibernetic -export_vtk
```

For each saved timestep number N a file `state_N.vtp` is created in the directory `./buffers`. Storing interval is given by the parameter `logstep`.

Run with Sibernetic-NEURON bridge
---------------------------------

Now it's possible to run the physical and neuronal simulations together. For this you need 
[sibernetic_NEURON](https://github.com/openworm/sibernetic_NEURON) also.  Don't forget to 
add the path of sibernetic_NEURON into your `PYTHONPATH`. You just need to run Sibernetic with 
command argument '-nrn <value>' where value is the path to NEURON simulation file (*.hoc e.g.). 
After that Sibernetic will initialise sibernetic_NEURON with the appropriate simulation file and 
same timeStep also. You should indicate from what segments of NEURON's model you'd like to read 
data (currently Voltage). After each step of the Sibernetic simulation it will run one step of 
the NEURON simulation and read data from it and update the signal array in Sibernetic. For now, 
it actually works in test mode list of segments is [hardcoded](https://github.com/openworm/sibernetic/blob/development/src/owNeuronSimulator.cpp#L70) 
so if you'd like to work with another list of segments you need rewrite this part of code and 
recompile Sibernetic.

If you have Sibernetic and [NEURON (with Python support)](http://neuralensemble.org/docs/PyNN/installation.html#installing-neuron) 
correctly installed, the following should be sufficient to get this running:

    git clone https://github.com/openworm/sibernetic_NEURON.git
    export PYTHONPATH=./sibernetic_NEURON:./src
    ./Release/Sibernetic -nrn ./sibernetic_NEURON/models/celegans/_ria.hoc  -f worm


Run with c302
---------------------------------

You can run Sibernetic with [c302](https://github.com/openworm/CElegansNeuroML/blob/master/CElegans/pythonScripts/c302/README.md) 
providing the input which will drive the contraction of the muscle cells.

If you have Sibernetic, [NEURON (with Python support)](http://neuralensemble.org/docs/PyNN/installation.html#installing-neuron) 
and [pyNeuroML](https://github.com/NeuroML/pyNeuroML) correctly installed, the following should be sufficient to get this running:

    git clone https://github.com/openworm/CElegansNeuroML.git
    export C302_HOME=./CElegansNeuroML/CElegans/pythonScripts/c302
    export PYTHONPATH=$PYTHONPATH:$C302_HOME:./src
    python sibernetic_c302.py

This will generate the NEURON code for the c302 simulation (using pyNeuroML), run Sibernetic with the neuronal simulation of c302 running in 
Python Neuron in the background, and save the results to files in the *simulations* directory (no Sibernetic gui will be shown). The 
simulation can be rerun with:

    ./Release/Sibernetic -l_from lpath=simulations/SimulationName_SimulationDate

For more information on options type:

    python sibernetic_c302.py -?



Making videos (*nix)
--------------------
If you run a simulation you may be interested in recording the graphical output.
You can either save the results to VTK files and use Paraview to create the video, or create the video using the default OpenGL visualisation.
Making such videos is a bit tricky because they need to be speeded up, so far I have found the following two commands do a decent job (change folder names accordingly) after you have used a screen record program:


If your video is in OGV format (if you used [recordmydesktop](http://recordmydesktop.sourceforge.net/about.php) for instance),
use the following script to convert to avi:

```
#!/bin/bash
 # ogv to avi
 # Call this with multiple arguments
 # for example : ls *.{ogv,OGV} | xargs ogv2avi
 N=$#;
 echo "Converting $N files !"
 for ((i=0; i<=(N-1); i++))
 do
 echo "converting" $1
 filename=${1%.*}
 mencoder "$1" -ovc xvid -oac mp3lame -xvidencopts pass=1 -o $filename.avi
 shift 1
 done
```

```
#make images from video
ffmpeg -i crawley_6.avi -r 0.05 -f image2 ~/Documents/tmp/output-%06d.jpg
```

```
#re-encode into video
ffmpeg -r 100 -i output-%06d.jpg -r 100 -vb 60M speeded_worm.mp4
```

Troubleshooting
--------------------
If you have any question or have a problem with running Sibernetic please contact with us.
Email me on skhayrulin@openworm.org or info@openworm.org. Or you can create an [issue on GitHub](https://github.com/openworm/sibernetic/issues).
