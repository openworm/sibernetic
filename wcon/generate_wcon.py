import math
import sys
import os
import numpy as np
import time


def dist(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def print_(msg):

    print("WCON gen >>> %s"%(msg))

def generate_wcon(pos_file_name,
                  wcon_file_name,
                  rate_to_plot=1,
                  max_lines=1e12,
                  max_time_s=1e12,
                  plot=False,
                  save_figure1_to=None,
                  save_figure2_to=None,
                  save_figure3_to=None):

    print_("Generating WCON from %s to %s, with plotting rate %i"%(pos_file_name,wcon_file_name,rate_to_plot))

    postions_file = open(pos_file_name)

    line_num = 0


    ts = []
    x = {}
    y = {}
    z = {}

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)

    num_plotted_frames = 0
    points_plotted = 0

    middle_points = []
    ave_points = []
    middle_point_speed_x = []
    middle_point_speed_y = []
    ave_point_speed_x = []
    ave_point_speed_y = []
    time_points = []

    num_frames = 0

    angles = {}

    xs0 = None
    ys0 = None
    for line in postions_file:

        line_num += 1
        if line_num > max_lines:
            print_("Finished parsing file, as max number of lines reached!")
            break

        if line_num == 1 or line_num % rate_to_plot == 0:
            num_frames += 1
            words = line.split()
            points = (len(words) - 4) / 3
            middle_point = points / 2
            t_s = float(words[0])
            ts.append(t_s)

            if t_s > max_time_s:
                print_("Finished parsing file, as max time reached!")
                break
            x[t_s] = [float(w) for w in words[2:2 + int(points)]]
            y[t_s] = [float(w) for w in words[3 + int(points):3 + 2 * int(points)]]
            z[t_s] = [float(w) for w in words[4 + 2 * int(points):]]

            print_("L%i: at time: %s sec found %i points: [(%s,%s,%s),...,(%s,%s,%s)]" % (
                line_num, t_s, points, x[t_s][0], y[t_s][0], z[t_s][0], x[t_s][int(points)-1], y[t_s][int(points)-1], z[t_s][int(points)-1]))

            spacing = 0.05

            xs = []
            ys = []
            avx = 0
            avy = 0
            offset = num_plotted_frames * spacing

            for i in range(len(x[t_s])):

                # Swap x and y so worm moves "up"
                xs.append(y[t_s][i] + offset)
                ys.append(x[t_s][i])

                avx += xs[-1] - offset
                avy += ys[-1]
                points_plotted += 1

            if xs0==None:
                xs0 = xs
            if ys0==None:
                ys0 = ys

            avx = avx / points
            avy = avy / points

            if len(middle_points) > 0:
                dt = t_s - time_points[-1]
                middle_point_speed_x.append(
                    (xs[int(middle_point)] - offset -
                     middle_points[-1][0]) / dt)

                middle_point_speed_y.append(
                    (ys[int(middle_point)] - middle_points[-1][1]) / dt)
                dav = dist(avx, avy, ave_points[-1][0], ave_points[-1][1])

                ave_point_speed_x.append((avx - ave_points[-1][0]) / dt)
                ave_point_speed_y.append((avy - ave_points[-1][1]) / dt)

                print_("  Speed of point" +
                      " %i: (%s,%s) -> (%s,%s): x %sum/s, y %sum/s" %
                      (middle_point, middle_points[-1][0],
                       middle_points[-1][1], xs[int(middle_point)],
                       ys[int(middle_point)], middle_point_speed_x[-1],
                       middle_point_speed_y[-1]))

                print_("  Speed of av point " +
                      "(%s,%s) -> : (%s,%s): x %sum/s, y %sum/s" %
                      (ave_points[-1][0], ave_points[-1][1], avx, avy,
                       ave_point_speed_x[-1], ave_point_speed_y[-1]))

            middle_points.append(
                (xs[int(middle_point)] - offset, ys[int(middle_point)]))
            ave_points.append((avx, avy))
            time_points.append(t_s)

            print_("  Plot frame %i at %s s; l %i: [(%s,%s),...#%i]" % (
                num_plotted_frames, t_s, line_num, xs[0], ys[0], len(xs)))

            num_plotted_frames += 1

            l = ax.plot(xs, ys, '-')
            if num_plotted_frames % 5 == 1:
                time_ = '%ss' % t_s if not t_s == int(t_s) \
                    else '%ss' % int(t_s)

                ax.text(50 + ((num_plotted_frames - 1) * spacing),
                        10, time_, fontsize=12)
            ax.plot([xx+offset for xx in xs0], ys0,':',color=l[0].get_color(),linewidth=0.5)

    data = np.zeros((len(ts),len(xs)-4))

    print_('Finished parsing %i lines'%line_num)

    for ti in range(len(ts)):

        for i in range(len(y[0]))[2:-2]:

            x1 = x[ts[ti]][i-2]
            xc = x[ts[ti]][i]
            x2 = x[ts[ti]][i+2]

            y1 = y[ts[ti]][i-2]
            yc = y[ts[ti]][i]
            y2 = y[ts[ti]][i+2]

            a1 = math.atan2(y1-yc,x1-xc)
            a2 = math.atan2(y2-yc,x2-xc)
            angle = a2-a1
            if angle <0: angle+=2*math.pi
            if angle >2*math.pi: angle= 2*math.pi-angle

            deg = 360*(angle/(2*math.pi))

            data[ti][i-2]= deg

            ###print("At t=%s, i=%s: angle from between (%s,%s) - (%s,%s) - (%s,%s) = %s, %sdeg"%(ts[ti], i, x1,y1,xc,yc,x2,y2,angle,deg))

    info = "Loaded: %s points from %s, saving %i frames" % (
        line_num, pos_file_name, num_frames)
    print_(info)

    wcon = open(wcon_file_name, 'w')

    wcon.write('''{
    "metadata":{
        "who":"sibernetic_c302",
        "timestamp":"%s",
        "protocol":"Generated by Sibernetic & c302!"
    },
    "units":{ "t":"s",
              "x":"micrometers",
              "y":"micrometers"},
    "comment":"Saved from Sibernetic data.",
    "note":"%s",
    "data":[\n''' % (time.strftime("%Y-%m-%dT%H:%M:%S+00:00", time.gmtime()), info))

    wcon.write('''            {"id":"wormTest",
             "t":[ ''')
    for t in ts:
        wcon.write('%s' % (t))
        if t != ts[-1]:
            wcon.write(', ')

    wcon.write('],\n')

    wcon.write('''             "x":[ ''')
    for t in ts:
        wcon.write('%s' % (x[t]))
        if t != ts[-1]:
            wcon.write(',\n')
    wcon.write('],\n')

    wcon.write('''             "y":[ ''')
    for t in ts:
        wcon.write('%s' % (y[t]))
        if t != ts[-1]:
            wcon.write(',\n')
    wcon.write(']\n')


    wcon.write('\n}\n')
    wcon.write('\n]\n')
    wcon.write('\n}\n')

    wcon.close()
    postions_file.close()

    info = "Midline of worm through time (to %s seconds)"%time_points[-1]

    plt.xlabel("x direction")
    plt.ylabel("y direction")
    fig.canvas.set_window_title(info)
    plt.title(info)

    if save_figure1_to:
        plt.savefig(save_figure1_to)

    fig = plt.figure()
    info = "Speed of worm in x (lateral) & y (along body) directions"
    fig.canvas.set_window_title(info)
    plt.title(info)
    plt.xlabel("Time (s)")
    plt.ylabel("Speed")

    plt.plot(time_points[1:], ave_point_speed_x, 'red',
             label='Speed in x of average of %i points' % points)
    plt.plot(time_points[1:], middle_point_speed_x, 'pink',
             label='Speed in x dir of point %i/%i' % (middle_point,
                                                      points), linestyle='--')
    plt.plot(time_points[1:], ave_point_speed_y, 'blue',
             label='Speed in y of average of %i points' % points)
    plt.plot(time_points[1:], middle_point_speed_y, 'cyan',
             label='Speed in y dir of point %i/%i' % (middle_point,
                                                      points), linestyle='--')



    plt.plot([0,time_points[-1]],[0,0], 'grey', linestyle=':')

    plt.legend(loc=2, fontsize = 'x-small')

    if save_figure2_to:
        plt.savefig(save_figure2_to,bbox_inches='tight')

    fig = plt.figure()

    plt.xlabel("Time (s)")
    plt.ylabel("Percentage along worm")

    plot0 = plt.imshow(data.transpose(), interpolation='nearest', aspect='auto')
    ax = plt.gca();

    info = "Propagation of curvature along body of worm (180=straight)"
    fig.canvas.set_window_title(info)
    plt.title(info)

    xt = ax.get_xticks()
    #print(xt)
    time_ticks = [time_points[int(ti)] if (ti>=0 and ti<len(time_points)) else 0 for ti in xt]
    #print(time_ticks)
    ax.set_xticklabels(time_ticks)

    fig.colorbar(plot0)

    if save_figure3_to:
        plt.savefig(save_figure3_to,bbox_inches='tight')

    if plot:
        plt.show()

    return x,y,z,ts


def validate(wcon_file):
    import json, jsonschema

    wcon_schema = "wcon_schema.json"

    if not os.path.isfile("wcon_schema.json"):
        print_("Cannot validate file: %s!! WCON schema %s not found!!"%(wcon_file, wcon_schema))
        return

    # The WCON schema
    with open(wcon_schema, "r") as wcon_schema_file:
        schema = json.loads(wcon_schema_file.read())

    # Our example WCON file
    with open(wcon_file, 'r') as infile:
        serialized_data = infile.read()

    # Load the whole JSON file into a nested dict.
    w = json.loads(serialized_data)

    # Validate the raw file against the WCON schema
    jsonschema.validate(w, schema)

    print_("File %s is valid WCON!!"%wcon_file)

def transform(i,max_rad=0.03):
    return 1+max_rad*math.sin((math.pi/4) * (i)/24.0)

def get_color_for_fract(fract):

    if fract<0: fract = 0
    if fract>1: fract = 1
    minCol = [.4,0,.9]
    maxCol = [1,0.5,0]
    return (minCol[0] + fract*(maxCol[0] - minCol[0]),\
            minCol[1] + fract*(maxCol[1] - minCol[1]),\
            minCol[2] + fract*(maxCol[2] - minCol[2]))

'''
def get_rainbow_color_for_volts(fract):

    if fract<0: fract = 0.0
    if fract>1: fract = 1.0

    hue = (270 * (1-fract))
    return "pigment { color CHSL2RGB(<%f,1,0.5>) } // v = %f, fract = %f"%( hue , v, fract)'''

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import matplotlib
    import math
    validate("test.wcon")

    test_sim = 'C0_Muscles_2018-01-19_15-59-59'
    #test_sim = 'Sibernetic_2018-01-17_21-10-11'
    #test_sim = 'C0_TargetMuscle_2018-01-23_14-19-54'

    '''
    test_sim = 'C2_FW_2018-01-27_14-52-26'
    test_sim="/home/padraig/git/OpenWorm/output/C2_FW_2018-01-28_13-35-36"
    pos_file_name = "%s/worm_motion_log.txt"%test_sim
    small_file = "%s/worm_motion_log.wcon"%test_sim
    x,y,z,ts = generate_wcon(pos_file_name, small_file, rate_to_plot=2, plot=False)

    sys.path.append("..")
    from plot_positions import plot_muscle_activity
    musc_act_file = "%s/muscles_activity_buffer.txt"%test_sim
    activations, times = plot_muscle_activity(musc_act_file,0.005,1000, show_plot=False)

    print("Plotting %s sets of %s (%s) activation values on %s times"%(len(activations['MDR']),len(activations['MDR'][0]), len(times),len(ts)))
    for ti in range(len(ts)):
        t = ts[ti]
        print("Check time index %i: %s ms"%(ti,t))


        fig = plt.figure()
        info = "Pos at %sms"%t
        fig.canvas.set_window_title(info)
        plt.title(info)

        mx = y[t]
        my = x[t]
        for li in range(len(mx)-1):
            print('------------')
            mi = 23-int(li/4.1)
            mt = ti*2
            scale = 1.5
            actDL = activations['MDL'][mi][mt]*scale
            colorDL = get_color_for_fract(actDL)
            actDR = activations['MDR'][mi][mt]*scale
            colorDR = get_color_for_fract(actDR)
            actVL = activations['MVL'][mi][mt]*scale
            colorVL = get_color_for_fract(actVL)
            actVR = activations['MVR'][mi][mt]*scale
            colorVR = get_color_for_fract(actVR)
            print('Midpoint %i mapped to muscle %s: %s - %s; %s - %s'%(li,mi, actDL, colorDL, actVL, colorVL))


            max_rad1 = 3.5
            max_rad2 = 2
            widthL = transform(li,max_rad=max_rad1)
            widthL1 = transform(li+1,max_rad=max_rad1)
            widthR = transform(li,max_rad=max_rad2)
            widthR1 = transform(li+1,max_rad=max_rad2)
            x0=mx[li]
            y0=my[li]
            x1=mx[li+1]
            y1=my[li+1]
            a=x1-x0
            b=y1-y0



            newx0 = x0+widthL*b
            newy0 = y0-widthL*a
            newx1 = x1+widthL1*b
            newy1 = y1-widthL1*a
            plt.plot([newx0,newx1],[newy0,newy1],color=colorDL,linewidth=0,marker='o',markersize=3,alpha=0.4, markeredgewidth=0)

            newx0 = x0+widthR*b
            newy0 = y0-widthR*a
            newx1 = x1+widthR1*b
            newy1 = y1-widthR1*a
            plt.plot([newx0,newx1],[newy0,newy1],color=colorDR,linewidth=0,marker='o',markersize=3,alpha=0.4, markeredgewidth=0)

            newx0 = x0-widthL*b
            newy0 = y0+widthL*a
            newx1 = x1-widthL1*b
            newy1 = y1+widthL1*a
            plt.plot([newx0,newx1],[newy0,newy1],color=colorVL,linewidth=0,marker='o',markersize=3,alpha=0.4, markeredgewidth=0)

            newx0 = x0-widthR*b
            newy0 = y0+widthR*a
            newx1 = x1-widthR1*b
            newy1 = y1+widthR1*a
            plt.plot([newx0,newx1],[newy0,newy1],color=colorVR,linewidth=0,marker='o',markersize=3,alpha=0.4, markeredgewidth=0)

        plt.plot(mx,my, 'green')

        bound = .7
        plt.xlim([mx[49]-bound*1.5,mx[49]+bound*1.5])
        plt.ylim([my[49]-bound,my[49]+bound])
        ax = plt.gca();
        #ax.set_facecolor('black')
        ax.set_aspect('equal')

        matplotlib.rc('axes',edgecolor='w')

        for child in ax.get_children():
            if isinstance(child, matplotlib.spines.Spine):
                child.set_color('w')
        ax.xaxis.grid(True,color='grey')
        ax.yaxis.grid(True,color='grey')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        frame = str(int(t*10000))
        while len(frame)<7:
            frame='0'+frame
        print("Writing frame %i: %s"%(ti,frame))
        plt.savefig("Calcium_%s.png"%frame,facecolor='black',transparent=True,bbox_inches='tight')

    exit()'''


    if len(sys.argv) >1:
        pos_file_name = sys.argv[1]
    else:
        pos_file_name = '../buffers/worm_motion_log.txt'

    if len(sys.argv) >2:
        rate_to_plot = int(sys.argv[2])
    else:
        rate_to_plot=1

    #generate_wcon(pos_file_name, "sibernetic_test_full.wcon", rate_to_plot=5)
    #generate_wcon(pos_file_name, "sibernetic_test_small.wcon",
    #              rate_to_plot=20, max_time_s=20)


    small_file = "small.wcon"
    generate_wcon(pos_file_name, small_file, rate_to_plot=rate_to_plot, plot=True)
    validate(small_file)
