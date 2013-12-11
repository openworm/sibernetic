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
