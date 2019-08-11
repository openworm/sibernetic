""" Simple configuration generator
"""
import enum
import json
import math


class ParticleType(enum.Enum):
    LIQUID = 1
    BOUND = 3

MASS = 20.00e-13

def make_particle(pos, vel=None, type=ParticleType.LIQUID, mass = MASS):
   particle_dict = {}
   particle_dict["position"] = pos
   particle_dict["type"] = type.value
   particle_dict["velocity"] = vel or [0.0, 0.0, 0.0, 1.0]
   particle_dict["mass"] = mass
   particle_dict["density"] = 1000.0
   return particle_dict

def gen(start, end, step):
    ret = start
    while ret < end:
        yield ret
        ret += step

def calc_delta(params):
    x = [1, 1, 0, -1, -1, -1, 0, 1, 1, 1,  0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 2, -2, 0, 0,  0,  0,  0, 0]
    y = [0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1,  0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 0, 2, -2, 0, 0,  0,  0]
    z = [0,  0,  0,  0,  0,  0,  0,  0,  1, 1, 1, 1, 1, 1,  1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 2, -2, 1, -1]
    sum1_x = 0.0
    sum1_y = 0.0
    sum1_z = 0.0
    sum1 = 0.0
    sum2 = 0.0
    v_x = 0.0
    v_y = 0.0
    v_z = 0.0
    dist = 0.0
    particleRadius = math.pow(params["mass"] / params["rho0"], 1.0 / 3.0)
    h_r_2 = 0.0
    for i in range(len(x)):
        v_x = x[i] * 0.8 * particleRadius
        v_y = y[i] * 0.8 * particleRadius
        v_z = z[i] * 0.8 * particleRadius
        dist = math.sqrt(v_x * v_x + v_y * v_y + v_z * v_z)
        if dist <= params["h"] * params["simulation_scale"]:
            h_r_2 = math.pow((params["h"] * params["simulation_scale"] - dist), 2)
            sum1_x += h_r_2 * v_x / dist
            sum1_y += h_r_2 * v_y / dist
            sum1_z += h_r_2 * v_z / dist

            sum2 += h_r_2 * h_r_2

    sum1 = sum1_x * sum1_x + sum1_y * sum1_y + sum1_z * sum1_z
    result = 1.0 / (params["beta"] * params["grad_wspiky_coefficient"] * params["grad_wspiky_coefficient"] * (sum1 + sum2))
    return result

def gen_param_map(out, h, x_dim, y_dim, z_dim, mass=MASS):
    out["parameters"] = {}
    out["parameters"]["particles"] = 30
    out["parameters"]["h"] = h
    out["parameters"]["x_max"] = h * x_dim
    out["parameters"]["x_min"] = 0
    out["parameters"]["y_max"] = h * y_dim
    out["parameters"]["y_min"] = 0
    out["parameters"]["z_max"] = h * z_dim
    out["parameters"]["z_min"] = 0
    out["parameters"]["mass"] = MASS
    out["parameters"]["rho0"] = 1000.0
    out["parameters"]["time_step"] = 4.0 * 5.0e-06
    out["parameters"]["simulation_scale"] = 0.0037 * math.pow(mass,1.0/3.0)/math.pow(0.00025, 1.0/3.0)

    out["parameters"]["beta"] = (out["parameters"]["time_step"] ** 2) * (out["parameters"]["mass"] **2) * 2 /(out["parameters"]["rho0"] ** 2);
    out["parameters"]["wpoly6_coefficient"] = 315.0 / ( 64.0 * math.pi * math.pow(h * out["parameters"]["simulation_scale"], 9.0 ) )
    out["parameters"]["grad_wspiky_coefficient"] = -45.0 / ( math.pi * math.pow( (h * out["parameters"]["simulation_scale"]), 6.0 ) )
    out["parameters"]["divgrad_wviscosity_coefficient"] = -out["parameters"]["grad_wspiky_coefficient"]
    out["parameters"]["mass_mult_wpoly6_coefficient"] = out["parameters"]["mass"] * out["parameters"]["wpoly6_coefficient"]
    out["parameters"]["mass_mult_grad_wspiky_coefficient"] = out["parameters"]["mass"] * out["parameters"]["grad_wspiky_coefficient"]
    out["parameters"]["mass_mult_divgrad_viscosity_coefficient"] = out["parameters"]["mass"]  * out["parameters"]["divgrad_wviscosity_coefficient"]
    out["parameters"]["surf_tens_coeff"] = out["parameters"]["mass_mult_wpoly6_coefficient"] * out["parameters"]["simulation_scale"]
    out["parameters"]["delta"] = calc_delta(out["parameters"])

    out["parameters"]["gravity_x"] =  0.0
    out["parameters"]["gravity_y"] =  -9.8
    out["parameters"]["gravity_z"] = 0.0
    out["parameters"]["mu"] = 0.1 * 0.00004

def gen_model(x_dim, y_dim, z_dim, file_name="tmp"):
    particle_count = 0
    h = 3.34
    r0 = h * 0.5

    out = {}
    gen_param_map(out, h, x_dim, y_dim, z_dim)
    out["model"] = []

    #draw_bounds(out, x_dim, y_dim, z_dim, h, r0)
    draw_bounds(out, x_dim, y_dim, z_dim, h, r0)
    for x in gen(2 * r0, h * x_dim // 2 - r0, r0):
        for y in gen(2 * r0, h * y_dim - r0, r0):
            for z in gen(2 * r0, h * z_dim - r0, r0):
                out["model"].append(make_particle(pos=[x, y, z, 1.0], type=ParticleType.LIQUID, mass=MASS))

    fp = open(file_name, 'w')
    json.dump(out, fp)
    print(len(out['model']), "particles generated")
    fp.close()
    old_gen(out)

def old_gen(model, file_name="old_tmp"):
    fp = open(file_name, 'w')
    fp.write("[simulation box]\n")
    fp.write(f'{model["parameters"]["x_min"]}\n')
    fp.write(f'{model["parameters"]["x_max"]}\n')
    fp.write(f'{model["parameters"]["y_min"]}\n')
    fp.write(f'{model["parameters"]["y_max"]}\n')
    fp.write(f'{model["parameters"]["z_min"]}\n')
    fp.write(f'{model["parameters"]["z_max"]}\n')

    fp.write("[position]\n")
    position_str = ""
    velocity_str = ""
    for p in model['model']:
        position_str += f"{p['position'][0]}\t{p['position'][1]}\t{p['position'][2]}\t{p['type'] + 0.1}\n"
        velocity_str += f"{p['velocity'][0]}\t{p['velocity'][1]}\t{p['velocity'][2]}\t{p['type'] + 0.1}\n"
    fp.write(position_str)
    fp.write("[velocity]\n")
    fp.write(velocity_str)
    fp.write("[connection]\n")
    fp.write("[membranes]\n")
    fp.write("[particleMemIndex]\n")
    fp.write("[end]\n")
    fp.close()

def draw_bounds(particles, x_dim, y_dim, z_dim, h, r0):
    nx = int( x_dim * h / r0 )  # Numbers of boundary particles on X-axis
    ny = int( y_dim * h / r0 )  # Numbers of boundary particles on Y-axis
    nz = int( z_dim * h / r0 )  # Numbers of boundary particles on Z-axis
    # 1 - top and bottom
    for ix in range(nx):
        for iy in range(ny):
            if ( ( ix == 0 ) or ( ix == nx - 1) ) or ( (iy == 0) or (iy == ny - 1 ) ) :
                if ( ( ix == 0 ) or ( ix == nx - 1 ) ) and ( ( iy == 0 ) or ( iy == ny - 1 ) ): #corners
                    x = ix * r0 + r0 / 2.0
                    y = iy * r0 + r0 / 2.0
                    z = 0.0 * r0 + r0 / 2.0
                    vel_x = ( 1.0 * float( ix == 0 ) - 1.0 * float( ix == nx - 1 ) ) / math.sqrt( 3.0 )
                    vel_y = ( 1.0 * float( iy == 0 ) - 1.0 * float( iy == ny - 1 ) ) / math.sqrt( 3.0 )
                    vel_z = 1.0 / math.sqrt( 3.0 )
                    particles['model'].append(
                        make_particle(
                            pos=[x, y, z, 1.0],
                            vel=[vel_x, vel_y, vel_z, 1.0],
                            type=ParticleType.BOUND,
                        )
                    )

                    x = ix * r0 + r0 / 2.0
                    y = iy * r0 + r0 / 2.0
                    z = (nz - 1.0) * r0 + r0 / 2.0
                    vel_x = ( 1.0 * float( ix == 0 ) - 1.0 * float( ix == nx - 1 ) ) / math.sqrt( 3.0 )
                    vel_y = ( 1.0 * float( iy == 0 ) - 1.0 * float( iy == ny - 1 ) ) / math.sqrt( 3.0 )
                    vel_z = -1.0 / math.sqrt( 3.0 )
                    particles['model'].append(
                        make_particle(
                            pos=[x, y, z, 1.0],
                            vel=[vel_x, vel_y, vel_z, 1.0],
                            type=ParticleType.BOUND,
                        )
                    )
                else: #edges
                    x = ix * r0 + r0 / 2.0
                    y = iy * r0 + r0 / 2.0
                    z = 0.0 *r0 + r0 / 2.0
                    vel_x = ( 1.0 * ( float( ix == 0 ) - float( ix == nx - 1 ) ) ) / math.sqrt( 2.0 )
                    vel_y = ( 1.0 * ( float( iy == 0 ) - float( iy == ny - 1 ) ) ) / math.sqrt( 2.0 )
                    vel_z = 1.0 / math.sqrt( 2.0 )
                    particles['model'].append(
                        make_particle(
                            pos=[x, y, z, 1.0],
                            vel=[vel_x, vel_y, vel_z, 1.0],
                            type=ParticleType.BOUND,
                        )
                    )
                    x = ix * r0 + r0 / 2.0
                    y = iy * r0 + r0 / 2.0
                    z = (nz - 1.0) * r0 + r0 / 2.0
                    vel_x = ( 1.0 * ( float( ix == 0 ) - float( ix == nx - 1 ) ) ) / math.sqrt( 2.0 )
                    vel_y = ( 1.0 * ( float( iy == 0 ) - float( iy == ny - 1 ) ) ) / math.sqrt( 2.0 )
                    vel_z = -1.0 / math.sqrt( 2.0 )
                    particles['model'].append(
                        make_particle(
                            pos=[x, y, z, 1.0],
                            vel=[vel_x, vel_y, vel_z, 1.0],
                            type=ParticleType.BOUND,
                        )
                    )
            else: #planes
                x = ix * r0 + r0 / 2.0
                y = iy * r0 + r0 / 2.0
                z = 0.0 * r0 + r0 / 2.0
                vel_x = 0.0
                vel_y = 0.0
                vel_z = 1.0
                particles['model'].append(
                    make_particle(
                        pos=[x, y, z, 1.0],
                        vel=[vel_x, vel_y, vel_z, 1.0],
                        type=ParticleType.BOUND,
                    )
                )

                x = ix * r0 + r0 / 2.0
                y = iy * r0 + r0 / 2.0
                z = (nz - 1.0) * r0 + r0 / 2.0
                vel_x = 0.0
                vel_y = 0.0
                vel_z = -1.0
                particles['model'].append(
                    make_particle(
                        pos=[x, y, z, 1.0],
                        vel=[vel_x, vel_y, vel_z, 1.0],
                        type=ParticleType.BOUND,
                    )
                )
    #2 - side walls OX-OZ and opposite
    for ix in range(nx):
        for iz in range(1,nz - 1):
            if (ix == 0) or (ix == nx - 1):
                x = ix * r0 + r0 / 2.0
                y = 0.0 * r0 + r0 / 2.0
                z = iz * r0 + r0 / 2.0
                vel_x = 0.0
                vel_y = 1.0 / math.sqrt(2.0)
                vel_z = 1.0 * ( float( iz == 0 ) - float(iz == nz - 1 ) ) / math.sqrt(2.0)
                particles['model'].append(
                    make_particle(
                        pos=[x, y, z, 1.0],
                        vel=[vel_x, vel_y, vel_z, 1.0],
                        type=ParticleType.BOUND,
                    )
                )
                x = ix * r0 + r0 / 2.0
                y = ( ny - 1 ) * r0 + r0 / 2.0
                z = iz * r0 + r0 / 2.0
                vel_x = 0.0
                vel_y = -1.0 / math.sqrt(2.0)
                vel_z = 1.0 * (float(iz == 0) - float(iz == nz - 1)) / math.sqrt(2.0)

                particles['model'].append(
                    make_particle(
                        pos=[x, y, z, 1.0],
                        vel=[vel_x, vel_y, vel_z, 1.0],
                        type=ParticleType.BOUND,
                    )
                )

            else:
                x = ix * r0 + r0 / 2.0
                y = 0.0 * r0 + r0 / 2.0
                z = iz * r0 + r0 / 2.0
                vel_x = 0.0
                vel_y = 1.0
                vel_z = 0.0
                particles['model'].append(
                    make_particle(
                        pos=[x, y, z, 1.0],
                        vel=[vel_x, vel_y, vel_z, 1.0],
                        type=ParticleType.BOUND,
                    )
                )
                x = ix * r0 + r0 / 2.0
                y = ( ny - 1 ) * r0 + r0 / 2.0
                z = iz * r0 + r0 / 2.0
                vel_x = 0.0
                vel_y = -1.0
                vel_z = 0.0
                particles['model'].append(
                    make_particle(
                        pos=[x, y, z, 1.0],
                        vel=[vel_x, vel_y, vel_z, 1.0],
                        type=ParticleType.BOUND,
                    )
                )

    #3 - side walls OY-OZ and opposite
    for iy in range(1,ny - 1):
        for iz in range(1,nz - 1):
            x = 0.0 * r0 + r0 / 2.0
            y = iy * r0 + r0 / 2.0
            z = iz * r0 + r0 / 2.0
            vel_x = 1.0
            vel_y = 0.0
            vel_z = 0.0
            particles['model'].append(
                make_particle(
                    pos=[x, y, z, 1.0],
                    vel=[vel_x, vel_y, vel_z, 1.0],
                    type=ParticleType.BOUND,
                )
            )

            x = (nx - 1) * r0 + r0 / 2.0
            y = iy * r0 + r0 / 2.0
            z = iz * r0 + r0 / 2.0
            vel_x = -1.0
            vel_y = 0.0
            vel_z = 0.0

            particles['model'].append(
                make_particle(
                    pos=[x, y, z, 1.0],
                    vel=[vel_x, vel_y, vel_z, 1.0],
                    type=ParticleType.BOUND,
                )
            )

def main():
    #gen_model(12, 6, 12)
    gen_model(48, 12, 24)

    #gen_model(96, 24, 48)
    #gen_model(96, 48, 96)
    #gen_model(192, 48, 192)
    #gen_model(6, 6, 6)


if __name__ == '__main__':
    main()
