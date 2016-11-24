import math
import numpy as np
import sys

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

def print_(msg):
    pre = "Python >> "
    print('%s %s'%(pre,msg.replace('\n','\n'+pre)))

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
        
        self.contraction_array[0][muscle_row_count - 2] = 0
        self.contraction_array[0][muscle_row_count - 1] = 0
        self.contraction_array[1][muscle_row_count - 2] = 0
        self.contraction_array[1][muscle_row_count - 1] = 0

        return list(np.concatenate([self.contraction_array[0],
                                    self.contraction_array[1],
                                    self.contraction_array[1],
                                    self.contraction_array[0]])) 
                                    
    def save_results(self):
        
        print_("> NOT saving results for MuscleSimulation")
        
        
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
            
        print_("Loaded a list of %i activity traces (values %s->%s) at %i time points from %s"%(len(self.values[0]), min_, max_, len(self.values), activity_file))
        
        max_ = 4e-7
        
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
        print_("Returning %i values at time: %f s, step: %i (index %i): [%f, %f, %f, ...]"%(len(v), t, self.step, index, v[0], v[1], v[2]))
        #print v
        self.step += 1
        return list(v)  
    
    
    
class C302NRNSimulation():

    max_ca = 4e-7
    max_ca_found = -1
    
    def __init__(self, tstop=100, dt=0.005, activity_file=None, verbose=True):
        
        #from LEMS_c302_C1_Full_nrn import NeuronSimulation
        from LEMS_c302_nrn import NeuronSimulation
        
        import neuron
        self.h = neuron.h
        
        self.verbose = verbose
        
        self.ns = NeuronSimulation(tstop, dt)
        print_("Initialised C302NRNSimulation of length %s ms and dt = %s ms..."%(tstop,dt))
        
        
    def save_results(self):
        
        print_("> Saving results at time: %s"%self.h.t)
        
        self.ns.save_results()
        
    def run(self, skip_to_time=-1):
        
        print_("> Current NEURON time: %s ms"%self.h.t)
        
        self.ns.advance()
        
        print_("< Current NEURON time: %s ms"%self.h.t)
                  
        values = list([self._scale(self.h.a_MDR01[0].soma.cai,print_it=True), \
                 self._scale(self.h.a_MVR01[0].soma.cai), \
                 self._scale(self.h.a_MVL01[0].soma.cai), \
                 self._scale(self.h.a_MDL01[0].soma.cai), \
                 self._scale(self.h.a_MDR02[0].soma.cai), \
                 self._scale(self.h.a_MVR02[0].soma.cai), \
                 self._scale(self.h.a_MVL02[0].soma.cai), \
                 self._scale(self.h.a_MDL02[0].soma.cai), \
                 self._scale(self.h.a_MDR03[0].soma.cai), \
                 self._scale(self.h.a_MVR03[0].soma.cai), \
                 self._scale(self.h.a_MVL03[0].soma.cai), \
                 self._scale(self.h.a_MDL03[0].soma.cai), \
                 self._scale(self.h.a_MDR04[0].soma.cai), \
                 self._scale(self.h.a_MVR04[0].soma.cai), \
                 self._scale(self.h.a_MVL04[0].soma.cai), \
                 self._scale(self.h.a_MDL04[0].soma.cai), \
                 self._scale(self.h.a_MDR05[0].soma.cai), \
                 self._scale(self.h.a_MVR05[0].soma.cai), \
                 self._scale(self.h.a_MVL05[0].soma.cai), \
                 self._scale(self.h.a_MDL05[0].soma.cai), \
                 self._scale(self.h.a_MDR06[0].soma.cai), \
                 self._scale(self.h.a_MVR06[0].soma.cai), \
                 self._scale(self.h.a_MVL06[0].soma.cai), \
                 self._scale(self.h.a_MDL06[0].soma.cai), \
                 self._scale(self.h.a_MDR07[0].soma.cai), \
                 self._scale(self.h.a_MVR07[0].soma.cai), \
                 self._scale(self.h.a_MVL07[0].soma.cai), \
                 self._scale(self.h.a_MDL07[0].soma.cai), \
                 self._scale(self.h.a_MDR08[0].soma.cai), \
                 self._scale(self.h.a_MVR08[0].soma.cai), \
                 self._scale(self.h.a_MVL08[0].soma.cai), \
                 self._scale(self.h.a_MDL08[0].soma.cai), \
                 self._scale(self.h.a_MDR09[0].soma.cai), \
                 self._scale(self.h.a_MVR09[0].soma.cai), \
                 self._scale(self.h.a_MVL09[0].soma.cai), \
                 self._scale(self.h.a_MDL09[0].soma.cai), \
                 self._scale(self.h.a_MDR10[0].soma.cai), \
                 self._scale(self.h.a_MVR10[0].soma.cai), \
                 self._scale(self.h.a_MVL10[0].soma.cai), \
                 self._scale(self.h.a_MDL10[0].soma.cai), \
                 self._scale(self.h.a_MDR11[0].soma.cai), \
                 self._scale(self.h.a_MVR11[0].soma.cai), \
                 self._scale(self.h.a_MVL11[0].soma.cai), \
                 self._scale(self.h.a_MDL11[0].soma.cai), \
                 self._scale(self.h.a_MDR12[0].soma.cai), \
                 self._scale(self.h.a_MVR12[0].soma.cai), \
                 self._scale(self.h.a_MVL12[0].soma.cai), \
                 self._scale(self.h.a_MDL12[0].soma.cai), \
                 self._scale(self.h.a_MDR13[0].soma.cai), \
                 self._scale(self.h.a_MVR13[0].soma.cai), \
                 self._scale(self.h.a_MVL13[0].soma.cai), \
                 self._scale(self.h.a_MDL13[0].soma.cai), \
                 self._scale(self.h.a_MDR14[0].soma.cai), \
                 self._scale(self.h.a_MVR14[0].soma.cai), \
                 self._scale(self.h.a_MVL14[0].soma.cai), \
                 self._scale(self.h.a_MDL14[0].soma.cai), \
                 self._scale(self.h.a_MDR15[0].soma.cai), \
                 self._scale(self.h.a_MVR15[0].soma.cai), \
                 self._scale(self.h.a_MVL15[0].soma.cai), \
                 self._scale(self.h.a_MDL15[0].soma.cai), \
                 self._scale(self.h.a_MDR16[0].soma.cai), \
                 self._scale(self.h.a_MVR16[0].soma.cai), \
                 self._scale(self.h.a_MVL16[0].soma.cai), \
                 self._scale(self.h.a_MDL16[0].soma.cai), \
                 self._scale(self.h.a_MDR17[0].soma.cai), \
                 self._scale(self.h.a_MVR17[0].soma.cai), \
                 self._scale(self.h.a_MVL17[0].soma.cai), \
                 self._scale(self.h.a_MDL17[0].soma.cai), \
                 self._scale(self.h.a_MDR18[0].soma.cai), \
                 self._scale(self.h.a_MVR18[0].soma.cai), \
                 self._scale(self.h.a_MVL18[0].soma.cai), \
                 self._scale(self.h.a_MDL18[0].soma.cai), \
                 self._scale(self.h.a_MDR19[0].soma.cai), \
                 self._scale(self.h.a_MVR19[0].soma.cai), \
                 self._scale(self.h.a_MVL19[0].soma.cai), \
                 self._scale(self.h.a_MDL19[0].soma.cai), \
                 self._scale(self.h.a_MDR20[0].soma.cai), \
                 self._scale(self.h.a_MVR20[0].soma.cai), \
                 self._scale(self.h.a_MVL20[0].soma.cai), \
                 self._scale(self.h.a_MDL20[0].soma.cai), \
                 self._scale(self.h.a_MDR21[0].soma.cai), \
                 self._scale(self.h.a_MVR21[0].soma.cai), \
                 self._scale(self.h.a_MVL21[0].soma.cai), \
                 self._scale(self.h.a_MDL21[0].soma.cai), \
                 self._scale(self.h.a_MDR22[0].soma.cai), \
                 self._scale(self.h.a_MVR22[0].soma.cai), \
                 self._scale(self.h.a_MVL22[0].soma.cai), \
                 self._scale(self.h.a_MDL22[0].soma.cai), \
                 self._scale(self.h.a_MDR23[0].soma.cai), \
                 self._scale(self.h.a_MVR23[0].soma.cai), \
                 self._scale(self.h.a_MVL23[0].soma.cai), \
                 self._scale(self.h.a_MDL23[0].soma.cai), \
                 self._scale(self.h.a_MDR24[0].soma.cai), \
                 self._scale(0), \
                 self._scale(self.h.a_MVL24[0].soma.cai), \
                 self._scale(self.h.a_MDL24[0].soma.cai)])
                 
        if self.verbose:
            print_("Returning %s values: %s"%(len(values),values))
        return values
        
    
    def _scale(self,ca,print_it=False):
        
        self.max_ca_found = max(ca,self.max_ca_found)
        scaled = min(1,(ca/self.max_ca))
        if print_it: 
            print_("- Scaling %s to %s (max found: %s)"%(ca,scaled,self.max_ca_found))
        return scaled
        


if __name__ == '__main__':
    
    print_("This script is used by the Sibernetic C++ application")
    print_("Running it directly in Python will only plot the waves being generated for sending to the muscle cells...")
    
    try_c302 = '-c302dat' in sys.argv
    try_c302_nrn = '-c302nrn' in sys.argv
    testnrn = '-testnrn' in sys.argv
    
    skip_to_time = 0
    
    max_time = 1.0 # s
    
    time_per_step = 0.001  #  s
    increment = time_per_step/default_time_per_step
    num_plots = 1
    
    ms = MuscleSimulation(increment=increment)
    
    if try_c302:
        #ms = C302Simulation('configuration/test/c302/c302_B_Muscles.muscles.activity.dat', scale_to_max=True)
        ms = C302Simulation('configuration/test/c302/c302_C1_Muscles.muscles.activity.dat', scale_to_max=True)
        #ms = C302Simulation('../../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/TestMuscles.activity.dat')
        #ms = C302Simulation('../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/c302_B_Oscillator.muscles.activity.dat')
        skip_to_time = 0.05
        max_time = 0.2
        
    elif try_c302_nrn or testnrn:
        dt = 0.1  # ms
        max_time = .5  # s
        maxt = max_time*1000

        time_per_step = dt/1000  #  s
        increment = time_per_step/default_time_per_step

        ms = C302NRNSimulation(tstop=maxt, dt=dt, verbose=False) 
        
    
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
    

    if testnrn:
        for step in range(num_steps):

            ms.run()
            
        ms.save_results()

        quit()
        
    
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
            print_("At step %s (%s s)"%(step, t))
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
    
    
