import sys 

'''
    Plot the positions from a saved Sibernetic position_buffer.txt file
'''
def plot_positions(pos_file_name, rate_to_plot = 100, save_figure=True, show_plot=True):
    
    postions_file = open(pos_file_name)
    
    rate_to_plot = max(1,int(rate_to_plot))

    index = 0

    xmin=-1
    xmax=-1
    ymin=-1
    ymax=-1
    zmin=-1
    zmax=-1

    import matplotlib.pyplot as plt

    print("Loading: %s"%pos_file_name)

    fig = plt.figure()

    ax = fig.add_subplot(111)
    #plt.xlim([-350, 450]) # based on standard configuration
    #plt.ylim([-50, 750])

    max_lines = 2e121
    max_time_ms = 10000

    #colors = {2.1:'blue', 2.2:'red',1.1:'black',3.0:'green'}
    points_plotted = 0

    a = b = c = dt = steps_per_frame = 0

    frame = 0
    in_frame = 0
    xs = []
    ys = []
    
    num_plotted_frames = 0

    for line in postions_file:

        if index>max_lines:
            print("Finished parsing file, as max number of lines reached!")
            break

        if index==0: xmin=float(line)
        elif index==1: xmax=float(line)
        elif index==2: ymin=float(line)
        elif index==3: ymax=float(line)
        elif index==4: zmin=float(line)
        elif index==5: zmax=float(line)
        elif index==6: a=float(line)
        elif index==7: b=float(line)
        elif index==8: c=float(line)
        elif index==9: dt=float(line)
        elif index==10: steps_per_frame=float(line)

        elif index==11:
            x = [xmin,xmin,xmax,xmax,xmin]
            y = [zmin,zmax,zmax,zmin,zmin]
            #plt.plot(x,y,'-', color='grey')

        elif  index>=12:

            t_ms = frame * steps_per_frame * dt * 1000

            if t_ms>max_time_ms:
                print("Finished parsing file, as max time (%s ms) reached!"%max_time_ms)

                break

            plot_frame = (frame%rate_to_plot == 0)

            if plot_frame:
                w = line.split()
                m = float(w[3])
                if m>2 and m<3:
                    x = float(w[0])
                    y = float(w[1])
                    z = float(w[2])
                    #plt.plot(x,z,'.',color=colors[m])
                    xs.append(x+num_plotted_frames*30)
                    ys.append(z)
                    points_plotted+=1

            in_frame+=1
            if in_frame == a+b+c-1:

                if plot_frame:
                    print(" >> Plotting frame %i at %s ms; line %i: %s...\n"%(num_plotted_frames,t_ms,index,line))
                    ax.plot(xs,ys,'.', markersize=1)
                    num_plotted_frames+=1
                    if num_plotted_frames%3 == 1:
                        time = '%sms'%t_ms if not t_ms==int(t_ms) else '%sms'%int(t_ms)
                        ax.text(50+((num_plotted_frames-1)*30), 510, time, fontsize=12)

                frame+=1 
                in_frame = 0
                print("New positions (#%i) at time %s ms; line %i: %s"%(frame,t_ms,index,line))
                xs = []
                ys = []


        index+=1

    print("Loaded: %s points from %s, showing %s points in %i plots"%(index,pos_file_name,points_plotted,num_plotted_frames))

    if save_figure:
        plt.savefig('%s.png'%pos_file_name,bbox_inches='tight')
  
    if show_plot:
        plt.show()


def plot_muscle_activity(muscle_file_name, dt, logstep, show_plot=True):
    muscle_file = open(muscle_file_name)

    count = 0
    a = []
    times = []
    muscle_names = []
    for m in ['MDR','MVR','MVL','MDL']:
        for i in range(24):
            muscle_names.append('%s%s'%(m,i))
            
    for line in muscle_file:

        acts = [float(w) for w in line.split()]
        a.append(acts)
        times.append(count*dt*logstep)

        print("Found %s activation values (%s,...,%s) at line %s"%(len(acts),acts[0],acts[-1],count))
        count+=1
        

    import matplotlib.pyplot as plt
    import numpy as np
    
    aa = np.array(a).transpose()
    fig, ax = plt.subplots()
    plot0 = ax.pcolormesh(aa)
    ax.set_xticks(np.arange(aa.shape[1]) + 0.5, minor=False)
    ax.set_xticklabels(times)
    
    ax.set_yticks(np.arange(aa.shape[0]) + 0.5, minor=False)
    ax.set_yticklabels(muscle_names)
    
    fig.colorbar(plot0)
    
    if show_plot:
        plt.show()

if __name__ == '__main__':
    
    plot_muscle_activity('simulations/C0_TargetMuscle_2018-01-16_17-03-35/muscles_activity_buffer.txt',0.005,1000)
    exit()

    if len(sys.argv) == 2:
        pos_file_name = sys.argv[1]
    else:
        pos_file_name = 'buffers/position_buffer.txt'


    plot_positions(pos_file_name, rate_to_plot = 50, show_plot=True, save_figure=True)
