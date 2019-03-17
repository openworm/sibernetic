""" Simple configuration generator
"""


def gen(start, end, step):
    ret = start
    while ret < end:
        yield ret
        ret += step

def gen_model(x_dim, y_dim, z_dim):
    particle_count = 0
    h = 3.34
    r0 = h * 0.5
    out_file = open("tmp", "w")
    out_file.write("parameters[\n")
    out_file.write("\tparticles:30\n")
    out_file.write(f"\tx_max: {h * x_dim}\n")
    out_file.write("\tx_min: 0\n")
    out_file.write(f"\ty_max: {h * y_dim}\n")
    out_file.write("\ty_min: 0\n")
    out_file.write(f"\tz_max: {h * z_dim}\n")
    out_file.write("\tz_min: 0\n")
    out_file.write("\tmass: 20.00e-13\n")
    out_file.write("\trho0: 1000.0\n")
    out_file.write("\ttime_step: 20.0e-06\n")
    out_file.write("]\n")
    out_file.write("model[\n")
    out_file.write("\tposition[\n")
    for x in gen(r0, h * x_dim, r0):
        for y in gen(r0, h * y_dim, r0):
            for z in gen(r0, h * z_dim, r0):
                particle_count += 1
                out_file.write("\t\t{0} {1} {2} 1\n".format(x, y, z))
    out_file.write("\t]\n")
    out_file.write("\tvelocity[\n")
    for _ in range(particle_count):
        out_file.write("\t\t0 0 0 1\n")
    out_file.write("\t]\n")
    out_file.write("]\n")
    out_file.close()
    print(particle_count, "particles generated")


def main():
    gen_model(6, 6, 4)


if __name__ == '__main__':
    main()
