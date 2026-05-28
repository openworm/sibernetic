#This script will run a number of tests using sibernetic_c302.py & check the 
# files produced
set -e

# No c302
python sibernetic_c302.py -test -noc302 -duration 0.1

python sibernetic_c302.py -test -noc302 -duration 0.054 -logstep 3

# c302
python sibernetic_c302.py -test  -duration 1.1  -c302params C1 

# c302 + half_resolution
python sibernetic_c302.py -test  -duration 1  -c302params C0 -configuration worm_alone_half_resolution 

# c302 + TestMuscle 
python sibernetic_c302.py -test -duration 20 -c302params C0 -reference TargetMuscle -configuration worm_alone_half_resolution -logstep 500


