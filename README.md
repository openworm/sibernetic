![Sibernetic](http://i.imgur.com/Hbsw6Zs.png)

Sibernetic is physical simulator of biomechanical matter (membranes, elastic matter, contractile matter) and environments (liquids, solids and elastic matter with variable physical properties) developed for simulations of C. elegans physical body dynamics within the [OpenWorm project](http://www.openworm.org) by Andrey Palyanov, Sergey Khayrulin and Mike Vella (development of a Python module for external muscle activating signals generation and input) as part of the [OpenWorm team](http://www.openworm.org/people.html). At its core, Sibernetic is built as an extension to Predictive-Corrective Incompressible Smoothed Particle Hydrodynamics (PCISPH). It is primarily written in  C++ and OpenCL, which makes possible to run simulations on CPUs or GPUs, and has 3D visualization support built on top of OpenGL. 

There is a separate effort lead by [Giovanni Idili](https://github.com/gidili) and [Sergey Khayrulin](https://github.com/skhayrulin) to port this code to Java, as part of the [Geppetto simulation framework](http://www.geppetto.org/). 

Compiling / running (Linux/mac)
------------------------------

[![Join the chat at https://gitter.im/openworm/sibernetic](https://badges.gitter.im/openworm/sibernetic.svg)](https://gitter.im/openworm/sibernetic?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

**Linux**

Install OpenCL on Ubuntu. We suggest you initially go with [AMD OpenCL drivers](http://developer.amd.com/tools-and-sdks/heterogeneous-computing/amd-accelerated-parallel-processing-app-sdk/downloads/) as we have found these to be the most stable and complete. You can also try [Intel's drivers](http://develnoter.blogspot.co.uk/2012/05/installing-opencl-in-ubuntu-1204.html). This step often causes problems, contact the [openworm-discuss](https://groups.google.com/forum/#!forum/openworm-discuss) mailing list if you encounter issues. The AMD drivers include samples in /opt/AMDAPP/samples/opencl/bin which you can use to verify your OpenCL support is working.

You'll also need a variety of libraries. In ubuntu, install the dependencies with:

```
sudo apt-get install g++ python-dev freeglut3-dev nvidia-opencl-dev libglu1-mesa-dev libglew-dev python-numpy
```

Next, navigate to the `Release` folder and run:

```
make clean
make all
```

**Mac**: stay in the top-level folder. You need before run export several environment variables:
 
```
export PYTHONHEADERDIR=/usr/local/Cellar/python/<version_of_installed_pythonFramework>/Python.framework/Headers/
export PYTHONLIBDIR=/usr/local/lib/python2...
export PYTHONFRAMEWORKDIR=/usr/local/Frameworks/
```
Then
```
make clean -f makefile.OSX
make all -f makefile.OSX
```

You should see an output which looks something like this:

```
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

```
export PYTHONPATH=$PYTHONPATH:'./src'
```

Finally, to run, run the command:

**Linux**:
```
./Release/Sibernetic
```
**Mac**:
```
./Release/Sibernetic
```

You may need to make `./Release/Sibernetic` executable like so:

```
chmod +x ./Release/Sibernetic
```

If you do not run from the top-level folder you will see an error which looks something like this:

```
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

There are two demo scenes generated for Sibernetic. The first one contains an elastic cube covered with liquid-impermeable membranes and liquid inside. The second one contains two elastic membranes attached to boundary (one of it has a liquid-impermeable membranes covering and another one hasn't such). 

The second one contains two elastic membranes attached to a boundary (one of them has liquid-impermeable membranes covering them and the other one doesn't).

To switch between demos you need to press the 1 or 2 keys respectively. To pause simulation you may press space bar.

References

1. B. Solenthaler, Predictive-Corrective Incompressible SPH. ACM Transactions on Graphics (Proceedings of SIGGRAPH), 28(3), 2009. 
2. M. Ihmsen, N. Akinci, M. Gissler, M. Teschner, Boundary Handling and Adaptive Time-stepping for PCISPH Proc. VRIPHYS, Copenhagen, Denmark, pp. 79-88, Nov 11-12, 2010.
3. M. Becker, M. Teschner. Weakly compressible SPH for free surface flows // Proceedings of the 2007 ACM SIGGRAPH/Eurographics symposium on Computer animation, pages 209-217.

Main command options
--------------
To start Sibernetic with argument print in command prompt next ./Release/Sibernetic -whatever
Available options:
```
 -g_no                 Run without graphics
 -l_to                 Save simulation results to disk.
     logstep=<value>   Log every <value> steps
 -l_from               Load simulation results from disk.
     lpath=<value>     Indicates path where result of simulation will be stored. 
                       This option work only for -l_to and -l_from options
 -test                 Run some physical tests.
 -f <filename>         Load configuration from file <filename>.
 device=<device_type>  Trying to init OpenCL on device <type> it could be cpu or gpu 
                       default-ALL (it will try to init most powerful available device).
 timestep=<value>      Start simulation with time step = <value> in seconds.
 timelimit=<value>     Run simulation until <value> will be reached in seconds.
 leapfrog              Use for integration LeapFrog method
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
```
First 6 lines is a spatial description of boundary box
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

This create 3 new files in the folder ./buffers:
- connection_buffers.txt - stores information about connection among the elastic particles
- membranes_buffer.txt   - stores information about membranes 
- position_buffer.txt    - stores information about current position of all of the non boundary particles it save information to this file every 10 steps of simulation. You should remember that the more info you 
want to store than bigger output file is. 

For view result you should run simulation with:
```
./Release/Sibernetic -l_from
```
It get positions from position_buffer.txt file and displays the evolution of system in time


Making videos (*nix)
--------------------
If you run a simulation you may be interested in recording the graphical output. Making such videos is a bit tricky 
because they need to be speeded up, so far I have found the following two commands do a decent job (change folder 
names accordingly) after you have used a screen record program:


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
