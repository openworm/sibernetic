import math
import numpy as np

import matplotlib.pyplot as plt
from pylab import *

muscle_row_count = 24

default_time_per_step = 0.000001  #  s  
time_per_step = default_time_per_step  #  s 

quadrant0 = 'MDR'
quadrant1 = 'MVR'
quadrant2 = 'MVL'
quadrant3 = 'MDL'

colours = {}
colours[quadrant0] = '#000000'
colours[quadrant1] = '#00ff00'
colours[quadrant2] = '#0000ff'
colours[quadrant3] = '#ff0000'

"""

Get list of muscle names in same order as waves generated below. 
Based on info here:
https://github.com/openworm/Smoothed-Particle-Hydrodynamics/blob/3da1edc3b018c2e5c7c1a25e2f8d44b54b1a1c47/src/owWorldSimulation.cpp#L475

"""
def get_muscle_names():
    names = []
    for i in range(muscle_row_count):
        names.append(get_muscle_name(quadrant0, i))
        names.append(get_muscle_name(quadrant1, i))
        names.append(get_muscle_name(quadrant2, i))
        names.append(get_muscle_name(quadrant3, i))
    
    return names

def get_muscle_name(quadrant, index):
    return "%s%s"%(quadrant, index+1 if index>8 else ("0%i"%(index+1)))

def parallel_waves(n=muscle_row_count, #24 for our first test?
                   step=0, 
                   phi=math.pi,
                   amplitude=1,
				   #velocity=0.000008):
                   velocity =0.000015):
    """
    Array of two travelling waves, second one starts
    half way through the array
    """

    if n % 2 != 0:
        raise NotImplementedError("Currently only supports even number of muscles!")

    j = n/2

    row_positions = np.linspace(0,0.75*1.5*2*math.pi,j)

    wave_1 = (map(math.sin,(row_positions - velocity*step)))
    wave_2 = (map(math.sin,(row_positions + (math.pi) - velocity*step)))

    normalize_sine = lambda x : (x + 1)/2
    wave_1 = map(normalize_sine, wave_1)
    wave_2 = map(normalize_sine, wave_2)

    double_wave_1 = []
    double_wave_2 = []

    for i in wave_1:
        double_wave_1.append(i)
        double_wave_1.append(i)

    for i in wave_2:
        double_wave_2.append(i)
        double_wave_2.append(i)
        
    return (double_wave_1,double_wave_2)

class MuscleSimulation():

    def __init__(self,increment=1.0):
        self.increment = increment
        self.step = 0


    def run(self, skip_to_time=0, do_plot = True):
        """
		if(iterationCount<400000)
		{

			//muscle_activation_signal_cpp[0*24+20] = 0;
			//muscle_activation_signal_cpp[0*24+21] = 0;
			muscle_activation_signal_cpp[0*24+22] = 0;
			muscle_activation_signal_cpp[0*24+23] = 0;

			//muscle_activation_signal_cpp[1*24+20] = 0;
			//muscle_activation_signal_cpp[1*24+21] = 0;
			muscle_activation_signal_cpp[1*24+22] = 0;
			muscle_activation_signal_cpp[1*24+23] = 0;

			//muscle_activation_signal_cpp[2*24+20] = 0;
			//muscle_activation_signal_cpp[2*24+21] = 0;
			muscle_activation_signal_cpp[2*24+22] = 0;
			muscle_activation_signal_cpp[2*24+23] = 0;

			//muscle_activation_signal_cpp[3*24+20] = 0;
			//muscle_activation_signal_cpp[3*24+21] = 0;
			muscle_activation_signal_cpp[3*24+22] = 0;
			muscle_activation_signal_cpp[3*24+23] = 0;
		}
        """
        self.contraction_array =  parallel_waves(step = self.step)
        self.step += self.increment
		# for reversal movment after 40000 steps it will switch sinusoid
        if (self.step>400000):
            self.increment = -1.0
        else:
            self.contraction_array[0][muscle_row_count - 2] = 0
            self.contraction_array[0][muscle_row_count - 1] = 0
            self.contraction_array[1][muscle_row_count - 2] = 0
            self.contraction_array[1][muscle_row_count - 1] = 0

        return list(np.concatenate([self.contraction_array[0],
                                    self.contraction_array[1],
                                    self.contraction_array[1],
                                    self.contraction_array[0]]))  
class C302Simulation():
    
    values = []

    def __init__(self, 
                 activity_file='configuration/test/c302/c302_C1_Muscles.muscles.activity.dat', 
                 dt=0.0001,
                 scale_to_max=True):
                     
        self.step = 0
        self.dt = dt
        data = open(activity_file, 'r')
        min_ = 1e19
        max_ = -1e19
        for line in data:
            vv = []
            vs = line.strip().split('\t')
            for v in vs:
                vf = float(v)
                vv.append(vf)
                min_ = min(min_,vf)
                max_ = max(max_,vf)
            self.values.append(vv)
            
        print("Loaded a list of %i activity traces (values %s->%s) at %i time points from %s"%(len(self.values[0]), min_, max_, len(self.values), activity_file))
        
        max_ = 8e-8
        
        if scale_to_max:
            vals_scaled = []
            for vv in self.values:
                vv2 = [v/max_ for v in vv] 
                
                vals_scaled.append(vv2)
                
            self.values = vals_scaled

    def run(self, skip_to_time=0.05):
        t = skip_to_time + self.step*time_per_step
        
        index = int(t/self.dt)
        
        if (index<len(self.values)):
            v = self.values[index][1:48]
            v.append(0)
            v.extend(self.values[index][48:])
        else:
            v = np.zeros(96)
        print("Returning %i values at time: %f s, step: %i (index %i): [%f, %f, %f, ...]"%(len(v), t, self.step, index, v[0], v[1], v[2]))
        #print v
        self.step += 1
        return list(v)  
        


if __name__ == '__main__':
    
    print("This script is used by the Sibernetic C++ application")
    print("Running it directly in Python will only plot the waves being generated for sending to the muscle cells...")
    
    skip_to_time = 0
    
    max_time = 0.4 # s
    max_time = 0.4 # s
    
    time_per_step = 0.001  #  s
    increment = time_per_step/default_time_per_step
    num_plots = 1
    
    try_c302 = True
    #try_c302 = False
    
    ms = MuscleSimulation(increment=increment)
    
    if try_c302:
        #ms = C302Simulation('configuration/test/c302/c302_B_Muscles.muscles.activity.dat', scale_to_max=True)
        ms = C302Simulation('configuration/test/c302/c302_C1_Muscles.muscles.activity.dat', scale_to_max=True)
        #ms = C302Simulation('../../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/TestMuscles.activity.dat')
        #ms = C302Simulation('../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/c302_B_Oscillator.muscles.activity.dat')
        skip_to_time = 0.05
        max_time = 0.2
        
    
    activation = {}
    row = '11'
    row_int=int(row)
    m0='%s%s'%(quadrant0,row)
    m1='%s%s'%(quadrant1,row)
    m2='%s%s'%(quadrant2,row)
    m3='%s%s'%(quadrant3,row)
    
    for m in get_muscle_names():
        activation[m] = []
    times = []
    
    num_steps = int(max_time/time_per_step)
    steps_between_plots = int(num_steps/num_plots)
    
    
    
    show_all = True
        
        
    for step in range(num_steps):
        t = step*time_per_step
        
        l = ms.run(skip_to_time=skip_to_time)
        
        for i in range(muscle_row_count):
            mq0=get_muscle_name(quadrant0, i)
            activation[mq0].append(l[i])
            mq1=get_muscle_name(quadrant1, i)
            activation[mq1].append(l[i+muscle_row_count])
            mq2=get_muscle_name(quadrant2, i)
            activation[mq2].append(l[i+muscle_row_count*2])
            mq3=get_muscle_name(quadrant3, i)
            activation[mq3].append(l[i+muscle_row_count*3])
            
        times.append(t)
        
        if step==0 or step%steps_between_plots == 0:
            print "At step %s (%s s)"%(step, t)
            if show_all:
                figV = plt.figure()
                figV.suptitle("Muscle activation waves at step %s (%s s)"%(step, t))
                plV = figV.add_subplot(111, autoscale_on=True)

                plV.plot(l[0:muscle_row_count], label='%s*'%quadrant0,  color=colours[quadrant0], linestyle='-', marker='o')
                plV.plot(l[muscle_row_count:2*muscle_row_count], label='%s*'%quadrant1,  color=colours[quadrant1], linestyle='-', marker='o')
                plV.plot(l[2*muscle_row_count:3*muscle_row_count], label='%s*'%quadrant2,  color=colours[quadrant2], linestyle='-')
                plV.plot(l[3*muscle_row_count:4*muscle_row_count], label='%s*'%quadrant3,  color=colours[quadrant3], linestyle='-')

                plV.legend()
    
    
    if show_all:
        fig0 = plt.figure()
        fig0.suptitle("Muscle activation waves of [%s, %s, %s, %s] vs time"%(m0,m1,m2,m3))
        pl0 = fig0.add_subplot(111, autoscale_on=True)
        pl0.plot(times, activation[m0], label=m0,  color=colours[quadrant0], linestyle='-')
        pl0.plot(times, activation[m1], label=m1, color=colours[quadrant1], linestyle='-')
        pl0.plot(times, activation[m2], label=m2, color=colours[quadrant2], linestyle='--')
        pl0.plot(times, activation[m3], label=m3, color=colours[quadrant3], linestyle='--')

        pl0.legend()
        
    
    f, a = plt.subplots(4, sharex=True, sharey=True)
        
    a[0].set_title(quadrant0)
    a[0].set_ylabel('muscle #')
    a[1].set_title(quadrant3)
    a[1].set_ylabel('muscle #')
    a[2].set_title(quadrant1)
    a[2].set_ylabel('muscle #')
    a[3].set_title(quadrant2)
    a[3].set_ylabel('muscle #')
    a[3].set_xlabel('time step (%s s -> %s s)'%(skip_to_time, max_time))
    
    arr0 = []
    arr1 = []
    arr2 = []
    arr3 = []
    for i in range(muscle_row_count):
        arr0.append(activation[get_muscle_name(quadrant0, i)])
        arr1.append(activation[get_muscle_name(quadrant3, i)])
        arr2.append(activation[get_muscle_name(quadrant1, i)])
        arr3.append(activation[get_muscle_name(quadrant2, i)])
        
    a[0].imshow(arr0, interpolation='none', aspect='auto')
    a[1].imshow(arr1, interpolation='none', aspect='auto')
    a[2].imshow(arr2, interpolation='none', aspect='auto')
    a[3].imshow(arr3, interpolation='none', aspect='auto')
    
    #plt.colorbar()
        

    plt.show()
    
    
