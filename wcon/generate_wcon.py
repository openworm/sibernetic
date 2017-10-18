import math
import sys
import numpy as np


def dist(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def generate_wcon(pos_file_name,
                  wcon_file_name,
                  rate_to_plot=1,
                  max_lines=1e12,
                  max_time_s=1e12,
                  plot=False):

    postions_file = open(pos_file_name)

    line_num = 0

    print("Loading: %s" % pos_file_name)

    ts = []
    x = {}
    y = {}
    z = {}

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)

    num_plotted_frames = 0
    points_plotted = 0

    all_xs = {}
    all_ys = {}
    middle_points = []
    ave_points = []
    middle_point_speed_x = []
    middle_point_speed_y = []
    ave_point_speed_x = []
    ave_point_speed_y = []
    time_points = []

    num_frames = 0
    
    angles = {}

    for line in postions_file:

        line_num += 1
        if line_num > max_lines:
            print("Finished parsing file, as max number of lines reached!")
            break

        if line_num == 1 or line_num % rate_to_plot == 0:
            num_frames += 1
            words = line.split()
            points = (len(words) - 4) / 3
            middle_point = points / 2
            t_s = float(words[0])
            ts.append(t_s)

            if t_s > max_time_s:
                print("Finished parsing file, as max time reached!")
                break
            x[t_s] = [float(w) for w in words[2:2 + points]]
            y[t_s] = [float(w) for w in words[3 + points:3 + 2 * points]]
            z[t_s] = [float(w) for w in words[4 + 2 * points:]]

            print("L%i: at time: %s sec found %i points: [(%s,%s,%s),...,(%s,%s,%s)]" % (
                line_num, t_s, points, x[t_s][0], y[t_s][0], z[t_s][0], x[t_s][points-1], y[t_s][points-1], z[t_s][points-1]))

            if plot:
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
                    
                avx = avx / points
                avy = avy / points

                if len(middle_points) > 0:
                    dt = t_s - time_points[-1]
                    middle_point_speed_x.append(
                        (xs[middle_point] - offset -
                         middle_points[-1][0]) / dt)

                    middle_point_speed_y.append(
                        (ys[middle_point] - middle_points[-1][1]) / dt)
                    dav = dist(avx, avy, ave_points[-1][0], ave_points[-1][1])

                    ave_point_speed_x.append((avx - ave_points[-1][0]) / dt)
                    ave_point_speed_y.append((avy - ave_points[-1][1]) / dt)

                    print("  Speed of point" +
                          " %i: (%s,%s) -> (%s,%s): x %sum/s, y %sum/s" %
                          (middle_point, middle_points[-1][0],
                           middle_points[-1][1], xs[middle_point],
                           ys[middle_point], middle_point_speed_x[-1],
                           middle_point_speed_y[-1]))

                    print("  Speed of av point " +
                          "(%s,%s) -> : (%s,%s): x %sum/s, y %sum/s" %
                          (ave_points[-1][0], ave_points[-1][1], avx, avy,
                           ave_point_speed_x[-1], ave_point_speed_y[-1]))

                middle_points.append(
                    (xs[middle_point] - offset, ys[middle_point]))
                ave_points.append((avx, avy))
                time_points.append(t_s)

                print("  Plot frame %i at %s s; l %i: [(%s,%s),...#%i]\n" % (
                    num_plotted_frames, t_s, line_num, xs[0], ys[0], len(xs)))
                    
                num_plotted_frames += 1
                
                ax.plot(xs, ys, '-')
                if num_plotted_frames % 5 == 1:
                    time = '%ss' % t_s if not t_s == int(t_s) \
                        else '%ss' % int(t_s)
    
                    ax.text(50 + ((num_plotted_frames - 1) * spacing),
                            10, time, fontsize=12)
                            
    data = np.zeros((len(ts),len(xs)-2))
    for ti in range(len(ts)):
        
        for i in range(len(y[0]))[1:-1]:
            
            x1 = x[ts[ti]][i-1]
            xc = x[ts[ti]][i]
            x2 = x[ts[ti]][i+1]
            
            y1 = y[ts[ti]][i-1]
            yc = y[ts[ti]][i]
            y2 = y[ts[ti]][i+1]
            
            a1 = math.atan2(y1-yc,x1-xc) 
            a2 = math.atan2(y2-yc,x2-xc)
            angle = a2-a1
            if angle <0: angle+=2*math.pi
            if angle >2*math.pi: angle= 2*math.pi-angle
            
            deg = 360*(angle/(2*math.pi))
            
            data[ti][i-1]= deg
            
            print("At t=%s, i=%s: angle from between (%s,%s) - (%s,%s) - (%s,%s) = %s, %sdeg"%(ts[ti], i, x1,y1,xc,yc,x2,y2,angle,deg))

    info = "Loaded: %s points from %s, saving %i frames" % (
        line_num, pos_file_name, num_frames)
    print(info)

    wcon = open(wcon_file_name, 'w')

    wcon.write('''{
    "metadata":{
        "who":"Open Worm",
        "timestamp":"2016-01-22T17:44:48",
        "protocol":"Testing Sibernetic!"
    },
    "units":{ "t":"s",
              "x":"um",
              "y":"um"},
    "comment":"Saved from Sibernetic data.",
    "note":"%s",
    "data":[\n''' % info)

    for t in ts:
        wcon.write('''
           {"id":"worm1", "t":[%s],
            "x":[%s],
            "y":[%s]
           }''' % (t, x[t], y[t]))
        if t == ts[-1]:
            wcon.write('\n        ]')
        else:
            wcon.write(',\n')

    wcon.write('\n}\n')
    wcon.close()
    postions_file.close()

    if plot:

        fig = plt.figure()
        plt.plot(time_points[1:], middle_point_speed_x, 'cyan',
                 label='Speed in x dir of point %i/%i' % (middle_point,
                                                          points))
        plt.plot(time_points[1:], middle_point_speed_y, 'red',
                 label='Speed in y dir of point %i/%i' % (middle_point,
                                                          points))
        plt.plot(time_points[1:], ave_point_speed_x, 'blue',
                 label='Speed in x of average of %i points' % points)
        plt.plot(time_points[1:], ave_point_speed_y, 'green',
                 label='Speed in y of average of %i points' % points)
        plt.legend()
        
        fig = plt.figure()
        plot0 = plt.imshow(data)
        
        fig.colorbar(plot0)
        plt.show()
        

def validate(wcon_file):
    import json, jsonschema

    # The WCON schema
    with open("wcon_schema.json", "r") as wcon_schema_file:
        schema = json.loads(wcon_schema_file.read())

    # Our example WCON file
    with open(wcon_file, 'r') as infile:
        serialized_data = infile.read()

    # Load the whole JSON file into a nested dict.
    w = json.loads(serialized_data)

    # Validate the raw file against the WCON schema
    jsonschema.validate(w, schema)
    
    print("File %s is valid WCON!!"%wcon_file)


if __name__ == '__main__':

    validate("test.wcon")
    if len(sys.argv) == 2:
        pos_file_name = sys.argv[1]
    else:
        pos_file_name = '../buffers/worm_motion_log.txt'

    #generate_wcon(pos_file_name, "sibernetic_test_full.wcon", rate_to_plot=5)
    #generate_wcon(pos_file_name, "sibernetic_test_small.wcon",
    #              rate_to_plot=20, max_time_s=20)

    
    small_file = "small.wcon"
    generate_wcon(pos_file_name, small_file, rate_to_plot=10, plot=True)
    validate(small_file)
    
    
