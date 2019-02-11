""" Simple configuration generator
"""


def gen(start, end, step):
    ret = start
    while ret < end:
        yield ret
        ret += step


def gen_model():
    particle_count = 0
    h = 3.34
    r0 = h * 0.5
    out_file = open("tmp", "w")
    out_file.write("model[\n")
    out_file.write("\tposition[\n")
    for x in gen(r0, h * 4, r0):
        for y in gen(r0, h * 4, r0):
            for z in gen(r0, h * 9, r0):
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
    gen_model()


if __name__ == '__main__':
    main()
