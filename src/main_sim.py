import math
import numpy as np

import matplotlib.pyplot as plt
from pylab import *

muscle_row_count = 24

time_per_step = 0.000005  #  ms

def parallel_waves(n=muscle_row_count, #26 for our first test?
                   step=0, 
                   phi=math.pi,
                   amplitude=1,
                   velocity=0.0001):
    """
    Array of two travelling waves, second one starts
    half way through the array
    """

    if n % 2 != 0:
        raise NotImplementedError("Currently only supports even number of muscles!")

    j = n/2

    row_positions = np.linspace(0,1.5*2*math.pi,j)

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

class muscle_simulation():

    def __init__(self,increment=1.0):
        self.increment = increment
        self.step = 0

    def run(self,do_plot = True):
        self.contraction_array =  parallel_waves(step = self.step)
        self.step += self.increment
        return list(np.concatenate([self.contraction_array[0],
                                    self.contraction_array[1],
                                    self.contraction_array[1],
                                    self.contraction_array[0]]))  
        #return(self.contraction_array)
        



if __name__ == '__main__':
    
    print("This script is used by the Sibernetic C++ application")
    print("Running it directly in Python will only plot the waves being generated for sending to the muscle cells...")
    
    ms = muscle_simulation()
    
    num_plots = 5
    steps = 20000
    m0 = []
    m1 = []
    times = []
    
    
    for step in range(num_plots*steps):
        t = step*time_per_step
        l = ms.run()
        m0.append(l[0])
        m1.append(l[muscle_row_count])
        times.append(t)
        if step==0 or step%steps == 0:
            print "At step %s (%s ms)"%(step, t)
            figV = plt.figure()
            figV.suptitle("Muscle activation waves at step %s (%s ms)"%(step, t))
            plV = figV.add_subplot(111, autoscale_on=True)
            plV.plot(l[0:muscle_row_count], solid_joinstyle ='round', solid_capstyle ='round', color='#ff0000', linestyle='-', marker='o')
            plV.plot(l[muscle_row_count:muscle_row_count*2], solid_joinstyle ='round', solid_capstyle ='round', color='#00ff00', linestyle='-', marker='o')
    
    
    fig0 = plt.figure()
    fig0.suptitle("Muscle activation waves vs time")
    pl0 = fig0.add_subplot(111, autoscale_on=True)
    pl0.plot(times, m0,solid_joinstyle ='round', solid_capstyle ='round', color='#ff0000', linestyle='-')
    pl0.plot(times, m1,solid_joinstyle ='round', solid_capstyle ='round', color='#ffff00', linestyle='-')
    #print m0
    #pl0.plot(m1.values(), m1.keys(), solid_joinstyle ='round', solid_capstyle ='round', color='#ffff00', linestyle='-')

    plt.show()
