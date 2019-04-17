""" Simple configuration generator
"""
import enum
import json
import math


class ParticleType(enum.Enum):
    LIQUID = 1
    BOUND = 3


def make_particle(pos, vel=None, type=ParticleType.LIQUID, mass = 1.0):
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

    #draw_bounds(out, x_dim, y_dim, z_dim, h, r0)
    draw_bounds(out, x_dim, y_dim, z_dim, h, r0)
    for x in gen(2 * r0, h * x_dim - r0, r0):
        for y in gen(2 * r0, h * y_dim - r0, r0):
            for z in gen(2 * r0, h * z_dim - r0, r0):
                out["model"].append(make_particle(pos=[x, y, z, 1.0], type=ParticleType.LIQUID))

    fp = open(file_name, 'w')
    json.dump(out, fp)
    print(len(out['model']), "particles generated")
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
    gen_model(12, 6, 12)


if __name__ == '__main__':
    main()
