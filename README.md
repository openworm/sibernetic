Smoothed-Particle-Hydrodynamics
===============================

This is a C++ implementation of the Contractile SPH (Electrofluid) algorithm applied to C. elegans locomotion. Electrofluid is a modification of the [Smoothed Particle Hydrodynamics](http://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) developed for the [OpenWorm](http://openworm.org) project by Andrey Palyanov, Sergey Khayrulin and Mike Vella.

There is a separate effort lead by [Giovanni Idili](https://github.com/gidili) to port this code to Java, as part of the [Geppetto simulation framework](https://github.com/openworm/OpenWorm/wiki/Geppetto--Overview). 

Compiling / running (Linux/mac)
------------------------------

**Linux**: navigate to the `Release` folder and run:

```
make clean
make all
```

**Mac**: stay in the top-level folder and run:

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

Linux:
```
./Release/Smoothed-Particle-Hydrodynamics
```
Mac:
```
./build/Smoothed-Particle-Hydrodynamics
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

Saving to disk
--------------

You may wish to save simulations to disk rather than visualise them (**WARNING**: This is buggy)

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
If you run a simulation you may be interested in recording the graphical output. Making such videos is a bit tricky because they need to be speeded up, so far I have found the following two commands do a decent job (change folder names accordingly) after you have used a screen record program:

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
ffmpeg -r 100 -i output-%06d.jpg -r 100 -vb 60M speeded_worm.mp4
```
