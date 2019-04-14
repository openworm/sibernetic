import json
import math

from enum import Enum

class LoadMode(Enum):
    const = 1
    pos = 2
    vel = 3
    other = 4

FILE_NAME = 'test_energy'

h = 3.34
r0 = h * 0.5

x_min = 0
x_max = 20 * h
y_min = 0
y_max = 30 * h
z_min = 0
z_max = 20 * h

def gen_param_map(out, h, mass=1.0):
    out["parameters"] = {}
    out["parameters"]["particles"] = 30
    out["parameters"]["x_max"] = x_max
    out["parameters"]["x_min"] = x_min
    out["parameters"]["y_max"] = y_max
    out["parameters"]["y_min"] = y_min
    out["parameters"]["z_max"] = z_max
    out["parameters"]["z_min"] = z_min
    out["parameters"]["mass"] = 20.00e-13
    out["parameters"]["rho0"] = 1000.0
    out["parameters"]["time_step"] = 0.0001
    out["parameters"]["simulation_scale"] = 0.0041 * math.pow(mass,1.0/3.0)/math.pow(0.00025, 1.0/3.0)
    out["parameters"]["surf_tens_coeff"] = 1
    out["parameters"]["delta"] = 1


def make_particle(pos, vel=None, type=1.0, mass = 1.0):
    particle_dict = {}
    particle_dict["position"] = pos
    particle_dict["type"] = type
    particle_dict["velocity"] = vel or [0.0, 0.0, 0.0, 1.0]
    particle_dict["mass"] = mass
    particle_dict["density"] = 1000.0
    return particle_dict


def main():
    f = open(FILE_NAME, 'r')
    const = {}
    position = []
    velocity = []
    mode = LoadMode.const
    currContainer = []
    for line in f:
        line = line.strip(' \n')
        if line == '[position]':
            currContainer = position
            continue
        elif line == '[velocity]':
            currContainer = velocity
            continue
        elif line == '[connection]':
            break
        currContainer.append(list(map(float, line.split('\t'))))
    f.close()
    particles = list(zip(position, velocity))

    out = {}
    print(len(particles))
    particle_count = 0


    out = {}
    gen_param_map(out, h)
    out["model"] = []
    for p in particles:
        out["model"].append(make_particle(pos=p[0], vel=p[1], type=int(p[0][3] + 0.1), mass=20.00e-13))

    fp = open('converted', 'w')
    json.dump(out, fp)
    print(particle_count, "particles generated")
    fp.close()


if __name__ == '__main__':
    main()