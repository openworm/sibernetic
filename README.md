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

Making videos
-------------

Making videos is a bit tricky because they need to be speeded up, so far I have found the following two commands do a decent job (change folder names accordingly) after you have used a screen record program:

```
#If your video is in OGV  format (if you used recordmydesktop for instance), use the following script to convert to avi:

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
ffmpeg -r 100 -i output-%06d.jpg -r 100 -pix_fmt yuv420p speeded_worm.mp4
```

Compiling/ running (Linux/mac)
------------------------------

Navigate to the `Release` folder and run:

```
make clean
make all
```

You should see an output which looks something like this:

```
Building file: ../src/PyramidalSimulation.cpp
Invoking: GCC C++ Compiler
g++ -I/usr/include/python2.7/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/PyramidalSimulation.d" -MT"src/PyramidalSimulation.d" -o "src/PyramidalSimulation.o" "../src/PyramidalSimulation.cpp"
../src/PyramidalSimulation.cpp: In member function ‘std::vector<float> PyramidalSimulation::run()’:
../src/PyramidalSimulation.cpp:73:54: warning: deprecated conversion from string constant to ‘char*’ [-Wwrite-strings]
Finished building: ../src/PyramidalSimulation.cpp
 
Building file: ../src/main.cpp
Invoking: GCC C++ Compiler
g++ -I/usr/include/python2.7/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/main.d" -MT"src/main.d" -o "src/main.o" "../src/main.cpp"
In file included from ../src/owPhysicsFluidSimulator.h:6:0,
                 from ../src/owWorldSimulation.h:16,
                 from ../src/main.cpp:1:
../src/owOpenCLSolver.h:4:0: warning: ignoring #pragma comment  [-Wunknown-pragmas]
../src/main.cpp:4:13: warning: ‘load_to_file’ initialised and declared ‘extern’ [enabled by default]
../src/main.cpp:5:13: warning: ‘load_from_file’ initialised and declared ‘extern’ [enabled by default]
Finished building: ../src/main.cpp
 
Building file: ../src/owHelper.cpp
Invoking: GCC C++ Compiler
g++ -I/usr/include/python2.7/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/owHelper.d" -MT"src/owHelper.d" -o "src/owHelper.o" "../src/owHelper.cpp"
../src/owHelper.cpp: In function ‘int generateInnerWormLiquid(int, int, float*, float*)’:
../src/owHelper.cpp:541:6: warning: unused variable ‘segmentsCount’ [-Wunused-variable]
../src/owHelper.cpp:543:8: warning: unused variable ‘coeff’ [-Wunused-variable]
../src/owHelper.cpp:550:8: warning: unused variable ‘value’ [-Wunused-variable]
../src/owHelper.cpp:556:8: warning: unused variable ‘beta’ [-Wunused-variable]
../src/owHelper.cpp: In static member function ‘static void owHelper::generateConfiguration(int, float*, float*, float*&, int*, int&, int&, int&, int&, int&, int*)’:
../src/owHelper.cpp:736:9: warning: unused variable ‘r2ij’ [-Wunused-variable]
../src/owHelper.cpp:737:9: warning: unused variable ‘dx2’ [-Wunused-variable]
../src/owHelper.cpp:737:13: warning: unused variable ‘dy2’ [-Wunused-variable]
../src/owHelper.cpp:737:17: warning: unused variable ‘dz2’ [-Wunused-variable]
../src/owHelper.cpp:967:8: warning: unused variable ‘k’ [-Wunused-variable]
../src/owHelper.cpp:953:7: warning: variable ‘array_j’ set but not used [-Wunused-but-set-variable]
../src/owHelper.cpp:958:9: warning: unused variable ‘m_number’ [-Wunused-variable]
../src/owHelper.cpp:702:8: warning: unused variable ‘x’ [-Wunused-variable]
../src/owHelper.cpp:702:10: warning: unused variable ‘y’ [-Wunused-variable]
../src/owHelper.cpp:702:12: warning: unused variable ‘z’ [-Wunused-variable]
../src/owHelper.cpp:712:6: warning: unused variable ‘nEx’ [-Wunused-variable]
../src/owHelper.cpp:713:6: warning: unused variable ‘nEy’ [-Wunused-variable]
../src/owHelper.cpp:714:6: warning: unused variable ‘nEz’ [-Wunused-variable]
../src/owHelper.cpp:715:6: warning: unused variable ‘nMuscles’ [-Wunused-variable]
../src/owHelper.cpp:716:6: warning: unused variable ‘nM’ [-Wunused-variable]
../src/owHelper.cpp:716:9: warning: unused variable ‘nMi’ [-Wunused-variable]
../src/owHelper.cpp:716:13: warning: unused variable ‘nMj’ [-Wunused-variable]
../src/owHelper.cpp:717:6: warning: variable ‘wormIndex_start’ set but not used [-Wunused-but-set-variable]
../src/owHelper.cpp:717:22: warning: variable ‘wormIndex_end’ set but not used [-Wunused-but-set-variable]
../src/owHelper.cpp: In static member function ‘static void owHelper::preLoadConfiguration()’:
../src/owHelper.cpp:1428:7: warning: unused variable ‘i’ [-Wunused-variable]
../src/owHelper.cpp: In static member function ‘static void owHelper::loadConfigurationFromFile(float*&, float*&, int*&, int)’:
../src/owHelper.cpp:1668:8: warning: unused variable ‘i’ [-Wunused-variable]
../src/owHelper.cpp: In function ‘int _Z17generateWormShelliiPfS_RiPi.constprop.41(int, float*, float*, int&, int*)’:
../src/owHelper.cpp:95:5: warning: ‘prevSlice_start’ may be used uninitialised in this function [-Wuninitialized]
../src/owHelper.cpp:487:9: warning: ‘prevSlice_pCount’ may be used uninitialised in this function [-Wuninitialized]
../src/owHelper.cpp: In function ‘int generateWormShell(int, int, float*, float*, int&, int*)’:
../src/owHelper.cpp:95:5: warning: ‘prevSlice_start’ may be used uninitialised in this function [-Wuninitialized]
../src/owHelper.cpp:487:9: warning: ‘prevSlice_pCount’ may be used uninitialised in this function [-Wuninitialized]
Finished building: ../src/owHelper.cpp
 
Building file: ../src/owOpenCLSolver.cpp
Invoking: GCC C++ Compiler
g++ -I/usr/include/python2.7/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/owOpenCLSolver.d" -MT"src/owOpenCLSolver.d" -o "src/owOpenCLSolver.o" "../src/owOpenCLSolver.cpp"
In file included from ../src/owOpenCLSolver.cpp:5:0:
../src/owOpenCLSolver.h:4:0: warning: ignoring #pragma comment  [-Wunknown-pragmas]
../src/owOpenCLSolver.cpp: In member function ‘void owOpenCLSolver::initializeOpenCL()’:
../src/owOpenCLSolver.cpp:104:17: warning: unused variable ‘clSelectedPlatformID’ [-Wunused-variable]
../src/owOpenCLSolver.cpp: In member function ‘unsigned int owOpenCLSolver::_runIndexPostPass()’:
../src/owOpenCLSolver.cpp:304:37: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
Finished building: ../src/owOpenCLSolver.cpp
 
Building file: ../src/owPhysicsFluidSimulator.cpp
Invoking: GCC C++ Compiler
g++ -I/usr/include/python2.7/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/owPhysicsFluidSimulator.d" -MT"src/owPhysicsFluidSimulator.d" -o "src/owPhysicsFluidSimulator.o" "../src/owPhysicsFluidSimulator.cpp"
In file included from ../src/owPhysicsFluidSimulator.h:6:0,
                 from ../src/owPhysicsFluidSimulator.cpp:1:
../src/owOpenCLSolver.h:4:0: warning: ignoring #pragma comment  [-Wunknown-pragmas]
In file included from /usr/include/python2.7/Python.h:8:0,
                 from ../src/PyramidalSimulation.h:3,
                 from ../src/owPhysicsFluidSimulator.cpp:5:
/usr/include/python2.7/pyconfig.h:1161:0: warning: "_POSIX_C_SOURCE" redefined [enabled by default]
/usr/include/features.h:164:0: note: this is the location of the previous definition
/usr/include/python2.7/pyconfig.h:1183:0: warning: "_XOPEN_SOURCE" redefined [enabled by default]
/usr/include/features.h:166:0: note: this is the location of the previous definition
../src/owPhysicsFluidSimulator.cpp:9:12: warning: ‘iterationCount’ initialised and declared ‘extern’ [enabled by default]
../src/owPhysicsFluidSimulator.cpp:10:12: warning: ‘numOfElasticConnections’ initialised and declared ‘extern’ [enabled by default]
../src/owPhysicsFluidSimulator.cpp:11:12: warning: ‘numOfLiquidP’ initialised and declared ‘extern’ [enabled by default]
../src/owPhysicsFluidSimulator.cpp:12:12: warning: ‘numOfElasticP’ initialised and declared ‘extern’ [enabled by default]
../src/owPhysicsFluidSimulator.cpp:14:12: warning: ‘numOfMembranes’ initialised and declared ‘extern’ [enabled by default]
../src/owPhysicsFluidSimulator.cpp: In member function ‘double owPhysicsFluidSimulator::simulationStep()’:
../src/owPhysicsFluidSimulator.cpp:135:58: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
Finished building: ../src/owPhysicsFluidSimulator.cpp
 
Building file: ../src/owWorldSimulation.cpp
Invoking: GCC C++ Compiler
g++ -I/usr/include/python2.7/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/owWorldSimulation.d" -MT"src/owWorldSimulation.d" -o "src/owWorldSimulation.o" "../src/owWorldSimulation.cpp"
In file included from ../src/owPhysicsFluidSimulator.h:6:0,
                 from ../src/owWorldSimulation.h:16,
                 from ../src/owWorldSimulation.cpp:1:
../src/owOpenCLSolver.h:4:0: warning: ignoring #pragma comment  [-Wunknown-pragmas]
../src/owWorldSimulation.cpp:202:3: warning: "/*" within comment [-Wcomment]
../src/owWorldSimulation.cpp: In function ‘void drawScene()’:
../src/owWorldSimulation.cpp:483:6: warning: unused variable ‘n’ [-Wunused-variable]
../src/owWorldSimulation.cpp: In function ‘void renderInfo(int, int)’:
../src/owWorldSimulation.cpp:646:9: warning: unused variable ‘temp_v’ [-Wunused-variable]
Finished building: ../src/owWorldSimulation.cpp
 
Building target: Smoothed-Particle-Hydrodynamics
Invoking: GCC C++ Linker
g++ -L/usr/lib -L/usr/lib/python2.7 -o "Smoothed-Particle-Hydrodynamics"  ./src/PyramidalSimulation.o ./src/main.o ./src/owHelper.o ./src/owOpenCLSolver.o ./src/owPhysicsFluidSimulator.o ./src/owWorldSimulation.o   -lOpenCL -lpython2.7 -lrt -lglut -lGL -lGLU
Finished building target: Smoothed-Particle-Hydrodynamics
```

Then navigate to the top-level folder in the hierarchy (e.g `Smoothed-Particle-Hydrodynamics`) and set your `PYTHONPATH`:

```
export PYTHONPATH=$PYTHONPATH:'./src'
```

Finally, to run, run the command:

```
./Release/Smoothed-Particle-Hydrodynamics
```

You may need to make `./Release/Smoothed-Particle-Hydrodynamics` executable like so:

```
chmod +x ./Release/Smoothed-Particle-Hydrodynamics
```

If you do not run from the top-level folder you will see an error which looks something like this:

```
Compilation failed: 
"/tmp/OCLQ1BaOw.cl", line 8: catastrophic error: cannot open source file
"src//owOpenCLConstant.h"
#include "src//owOpenCLConstant.h"
```



