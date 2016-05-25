This takes the output of a the C302 Muscle Program with C1 Parameters and translated it into something that Sibernetic can understand.

***How to make it work on your machine:***

1. install ipython - <code>sudo apt-get ipython</code>

2. install neuron for python, for this you will have to build from source: [instructions](http://www.tc.umn.edu/~haszx010/files/vpl_dbs_docs/Installation.html)

3. Make sure that you are able to use the command nrngui, which means that neuron is in your PATH, if it you can't add the bin folder in neuron to your path. e.g. export <code>PATH=$PATH:/path/to/neuron/bin</code>

4. If you get this error:
dlopen failed -
x86_64/.libs/libnrnmech.so: undefined symbol: nrnmpi_myid | Make sure that the neuron python libraries to you PYTHONPATH e.g. <code>export PYTHONPATH=$PYTHONPATH:/path/to/neuron/lib/python</code>

5. Make sure that in inc/owSignalSimulator simClassName = "C302RealTimeSimulation"

6. (optional) If you want to import test simulations, you need pynml: [instructions](https://github.com/NeuroML/pyNeuroML)


***If you want to test a new C302 input, based on calcium concentration:***

1. Go to the c302 library and create a new test folder that includes only the files you need for the new simulation. For example, the c302_C1_Muscles simulation required
LEMS_c302_C1_Muscles.xml c302_C1_Muscles.nml and cell_C.xml. (If you don't know which files you need the pynml compiler will tell you which ones you are missing).

2. In this test folder run the pynml command. e.g. <code>pynml LEMS_c302_C1_Muscles.xml -neuron</code>

3. Copy this folder into the sibernetic project

4. In the C302RealTime class in main_sim.py change the activity file to point the folder that contains the c302 simulation

5. Make sure you scale the calcium values by the proper factor so that it is between 0-1 for example, with the C1 parameters, I divide each number by 1.56005925429e-12. You
can use the get_max function in the C302NrnSimulator (make sure you only use it to get the value the first time as it is slow) or you can use a known value.


***Here is a link to a g2.2xlarge aws ec2 image:***

[instance](https://console.aws.amazon.com/ec2/v2/home?region=us-east-1#LaunchInstanceWizard:ami=ami-815dabec)

[spot request](https://console.aws.amazon.com/ec2/v2/home?region=us-east-1#LaunchInstanceWizard:ami=ami-815dabec;type=spot)

Here is how you can see the output on your computer:

1. On remotehost <code>sudo service lightdm stop</code>

2. On remotehost  <code>X &</code>

3. On remotehost <code>x11vnc -display :0 -forever &</code>

4. On localhost <code>vncviewer</code>, when it asks for server type in the public IP of the instance, this can be found on the EC2 Console.

5. The first time that you run <code>./Release/Sibernetic/</code> you must run with sudo. After the first time, do not run as sudo again.

6. You need to add an inbound traffic rule for your security group:

Type: Custom TCP Rule
Protocol: TCP
Port Range: 5900 - 5920
Source: Anywhere


WARNING: In ec2 you are charged a small amount for data transferred over the network. It ends up being a very small amount of money, but just be cognizant of how much the vncviewer is running i.e. if you are not running a visual application, you should close down vnc viewer

