import math
import numpy as np

"""
import matplotlib.pyplot as plt
from pylab import *
"""

muscle_row_count = 24

time_per_step = 0.000001  #  s

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
        names.append("%s%s"%(quadrant0, i+1 if i>8 else ("0%i"%(i+1))))
    for i in range(muscle_row_count):
        names.append("%s%s"%(quadrant1, i+1 if i>8 else ("0%i"%(i+1))))
    for i in range(muscle_row_count):
        names.append("%s%s"%(quadrant2, i+1 if i>8 else ("0%i"%(i+1))))
    for i in range(muscle_row_count):
        names.append("%s%s"%(quadrant3, i+1 if i>8 else ("0%i"%(i+1))))
    
    return names


def parallel_waves(n=muscle_row_count, #26 for our first test?
                   step=0, 
                   phi=math.pi,
                   amplitude=1,
				   #velocity=0.000008):
                   velocity_s =0.000015*3.7*1.94/1.76, #swimming
				   velocity_c =0.000015*0.72): #crawling // 0.9 // 0.65
	"""
	Array of two travelling waves, second one starts
    half way through the array
	"""
	j = n/2
	max_muscle_force_coeff = 1.0
	
	# "<" = first 6 seconds crawling, then swimming
	# ">" = first 6 seconds swimming, then crawling
	if (step>1200000): 
		velocity = 14 * 0.000015*0.72#crawling
		max_muscle_force_coeff = 1.0
		row_positions = np.linspace(0,2.97*math.pi,j)		
		wave_m = np.linspace(1,0.6,j) 		
	else:
		velocity = 14 * 0.000015*3.7#swimming
		max_muscle_force_coeff = 0.575
		row_positions = np.linspace(0,0.81*math.pi,j)
		#wave_m = [0.8,0.7,0.8,0.93,1.0,1.0,1.0,1.0,0.93,0.8,0.6,0.4] 		
		#wave_m = [0.8,0.7,0.8,0.93,1.0,1.0,1.0,1.0,0.93,0.8,0.65,0.5] 		
		#wave_m = [0.8,0.7,0.8,0.93,1.0,1.0,1.0,0.93,0.8,0.6,0.4,0.2] 		
		#wave_m = [0.81,0.90,0.97,1.00,0.99,0.95,0.88,0.78,0.65,0.50,0.33,0.15] 		
		#wave_m = np.linspace(1,1,j) 		
		wave_m = [0.81,0.90,0.97,1.00,0.99,0.95,0.88,0.78,0.65,0.53,0.40,0.25] #6
		#wave_m = [0.81,0.90,0.97,1.00,0.99,0.95,0.90,0.83,0.75,0.65,0.55,0.45] 	#7	

	if n % 2 != 0:
		raise NotImplementedError("Currently only supports even number of muscles!")

	wave_1 = (map(math.sin, (row_positions - velocity*step) ))
	wave_2 = (map(math.sin, (row_positions - velocity*step + (math.pi)) ))
	
	normalize_sine = lambda x : abs(x*(x>0))#(x + 1)/2
	wave_1 = map(normalize_sine, wave_1)
	wave_2 = map(normalize_sine, wave_2)

	###### sinusoidal signal correction ##################################
	#normalize_sine = lambda x : x*x
	#wave_1 = map(normalize_sine, wave_1)
	#wave_2 = map(normalize_sine, wave_2)

	wave_1 = map(lambda x,y: max_muscle_force_coeff*x*y, wave_1, wave_m)
	wave_2 = map(lambda x,y: max_muscle_force_coeff*x*y, wave_2, wave_m)
	###### smooth start###################################################
	if (step<(10000/14)):	
		normalize_sine = lambda x : x*step/(10000/14)
		wave_1 = map(normalize_sine, wave_1)
		wave_2 = map(normalize_sine, wave_2)	
	######################################################################
	
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

    def run(self,do_plot = True):
        self.contraction_array =  parallel_waves(step = self.step)
        self.step += self.increment
        #if (self.step>1000000):
        #    self.increment = -1.0
        #else:
		#    self.increment = self.increment

        return list(np.concatenate([self.contraction_array[0],
                                    self.contraction_array[1],
                                    self.contraction_array[1],
                                    self.contraction_array[0]]))  
class C302Simulation():
    
    values = []

    def __init__(self, activity_file='configuration/test/c302/c302_B_Muscles.muscles.activity.dat', dt=0.0001):
        self.step = 0
        self.dt = dt
        data = open(activity_file, 'r')
        for line in data:
            vv = []
            vs = line.strip().split('\t')
            for v in vs:
                vv.append(float(v))
            self.values.append(vv)
            
        print("Loaded a list of %i activity traces at %i time points"%(len(self.values[0]), len(self.values)))
        

    def run(self, skip_to_time=0):
        t = skip_to_time + self.step*time_per_step
        
        index = int(t/self.dt)
        
        if (index<len(self.values)):
            v = self.values[index][1:48]
            v.append(0)
            v.extend(self.values[index][48:])
        else:
            v = np.zeros(96)
        #print("Returning %i values at time: %f s, step: %i (index %i): [%f, %f, %f, ...]"%(len(v), t, self.step, index, v[0], v[1], v[2]))
        #print v
        self.step += 1
        return list(v)  
        


if __name__ == '__main__':
    
    print("This script is used by the Sibernetic C++ application")
    print("Running it directly in Python will only plot the waves being generated for sending to the muscle cells...")
    
    ms = MuscleSimulation()
    ms = C302Simulation('configuration/test/c302/c302_B_Muscles.muscles.activity.dat')
    #ms = C302Simulation('../../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/TestMuscles.activity.dat')
    #ms = C302Simulation('../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/c302_B_Muscles.muscles.activity.dat')
    
    max_time = 0.4 # s
    num_plots = 4
    
    activation = {}
    row = '02'
    row_int=int(row)
    m0='%s%s'%(quadrant0,row)
    m1='%s%s'%(quadrant1,row)
    m2='%s%s'%(quadrant2,row)
    m3='%s%s'%(quadrant3,row)
    activation[m0] = []
    activation[m1] = []
    activation[m2] = []
    activation[m3] = []
    times = []
    
    num_steps = int(max_time/time_per_step)
    steps_between_plots = int(num_steps/num_plots)
    
    
    l = ms.run(skip_to_time=0.03)
        
        
    for step in range(num_steps):
        t = step*time_per_step
        activation[m0].append(l[row_int])
        activation[m1].append(l[row_int+muscle_row_count])
        activation[m2].append(l[row_int+muscle_row_count*2])
        activation[m3].append(l[row_int+muscle_row_count*3])
        times.append(t)
        if step==0 or step%steps_between_plots == 0:
            print "At step %s (%s s)"%(step, t)
            figV = plt.figure()
            figV.suptitle("Muscle activation waves at step %s (%s s)"%(step, t))
            plV = figV.add_subplot(111, autoscale_on=True)
            
            plV.plot(l[0:muscle_row_count], label='%s*'%quadrant0,  color=colours[quadrant0], linestyle='-', marker='o')
            plV.plot(l[muscle_row_count:2*muscle_row_count], label='%s*'%quadrant1,  color=colours[quadrant1], linestyle='-', marker='o')
            plV.plot(l[2*muscle_row_count:3*muscle_row_count], label='%s*'%quadrant2,  color=colours[quadrant2], linestyle='-')
            plV.plot(l[3*muscle_row_count:4*muscle_row_count], label='%s*'%quadrant3,  color=colours[quadrant3], linestyle='-')
            
            plV.legend()
    
    
    fig0 = plt.figure()
    fig0.suptitle("Muscle activation waves vs time")
    pl0 = fig0.add_subplot(111, autoscale_on=True)
    pl0.plot(times, activation[m0], label=m0,  color=colours[quadrant0], linestyle='-')
    pl0.plot(times, activation[m1], label=m1, color=colours[quadrant1], linestyle='-')
    pl0.plot(times, activation[m2], label=m2, color=colours[quadrant2], linestyle='--')
    pl0.plot(times, activation[m3], label=m3, color=colours[quadrant3], linestyle='--')
    
    pl0.legend()

    plt.show()
    
    
