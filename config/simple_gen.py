""" Simple configuration generator
"""

import json
import math

def make_particle(pos, vel=None, type=1.0, mass = 1.0):
   particle_dict = {}
   particle_dict["position"] = pos
   particle_dict["type"] = type
   particle_dict["velocity"] = vel or [0.0, 0.0, 0.0, 1.0]
   particle_dict["mass"] = mass
   particle_dict["density"] = 1000.0
   return particle_dict

def gen(start, end, step):
    ret = start
    while ret < end:
        yield ret
        ret += step

def gen_param_map(out, h, x_dim, y_dim, z_dim, mass=1.0):
    out["parameters"] = {}
    out["parameters"]["particles"] = 30
    out["parameters"]["x_max"] = h * x_dim
    out["parameters"]["x_min"] = 0
    out["parameters"]["y_max"] = h * y_dim
    out["parameters"]["y_min"] = 0
    out["parameters"]["z_max"] = h * z_dim
    out["parameters"]["z_min"] = 0
    out["parameters"]["mass"] = 20.00e-13
    out["parameters"]["rho0"] = 1000.0
    out["parameters"]["time_step"] = 0.0001
    out["parameters"]["simulation_scale"] = 0.0041 * math.pow(mass,1.0/3.0)/math.pow(0.00025, 1.0/3.0)
    out["parameters"]["surf_tens_coeff"] = 1
    out["parameters"]["delta"] = 1

def gen_model(x_dim, y_dim, z_dim, file_name="tmp"):
    particle_count = 0
    h = 3.34
    r0 = h * 0.5

    out = {}
    gen_param_map(out, h, x_dim, y_dim, z_dim)
    out["model"] = []
    for x in gen(r0, h * x_dim, r0):
        for y in gen(r0, h * y_dim, r0):
            for z in gen(r0, h * z_dim, r0):
                particle_count += 1
                out["model"].append(make_particle(pos=[x, y, z, 1.0]))

    fp = open(file_name, 'w')
    json.dump(out, fp)
    print(particle_count, "particles generated")
    fp.close()


def main():
    gen_model(6, 6, 4)


if __name__ == '__main__':
    main()
