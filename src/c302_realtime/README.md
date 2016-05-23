This takes the output of a the C302 Muscle Program with C1 Parameters and translated it into something that Sibernetic can understand.

How it works:

Run pynml on LEMS_c302_C1_Muscles.xml -neuron in the c302 project, The output of this is put in configuration/tests/c302/muscles_c1


How to make it work on your machine:

1. install ipython (?) - sudo apt-get ipython
2. install neuron for python, for this you will have to build from source, [instructions](http://www.tc.umn.edu/~haszx010/files/vpl_dbs_docs/Installation.html)
3. Add to you PYTHONPATH your neuron python library eg: /usr/local/nrn/lib/python
