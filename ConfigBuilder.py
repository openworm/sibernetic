import os
import copy
import re

SIMULATION_BOX = "simulation box"
PHYSICAL_PARAMETERS = "physical parameters"
POSITION = "position"
VELOCITY = "velocity"
CONNECTION = "connection"
MEMBRANES = "membranes"
PARTICLE_MEM_INDEX = "particleMemIndex"

MAX_NEIGHBOR_COUNT = 32  # must match MAX_NEIGHBOR_COUNT in owOpenCLConstant.h

_PHYS_PARAM_RE = re.compile(r"^\s*(\w+)\s*:\s*([^\s/]+)")


class Configuration:
    """
    Class to hold configuration data
    """

    def __init__(self):
        self.simulation_box = []
        self.physical_parameters = {}
        self.particles = {}
        self.connections = []
        self.membranes = []
        self.particle_mem_index = []

    def is_liquid_type(self, particle_type):
        return int(particle_type) == 1

    def is_elastic_type(self, particle_type):
        return int(particle_type) == 2

    def is_boundary_type(self, particle_type):
        return int(particle_type) == 3

    def is_liquid(self, particle_id):
        if particle_id in self.particles:
            particle_type = self.particles[particle_id].get(POSITION)[3]
            return self.is_liquid_type(particle_type)
        return False

    def is_elastic(self, particle_id):
        if particle_id in self.particles:
            particle_type = self.particles[particle_id].get(POSITION)[3]
            return self.is_elastic_type(particle_type)
        return False

    def is_boundary(self, particle_id):
        if particle_id in self.particles:
            particle_type = self.particles[particle_id].get(POSITION)[3]
            return self.is_boundary_type(particle_type)
        return False

    def get_particles(
        self,
        translate=(0, 0, 0),
        include_liquid=True,
        include_elastic=True,
        include_boundary=True,
    ):
        particles = []

        for i in self.particles:
            if (
                (include_liquid and self.is_liquid(i))
                or (include_elastic and self.is_elastic(i))
                or (include_boundary and self.is_boundary(i))
            ):
                pos = self.particles[i].get(POSITION)
                vel = self.particles[i].get(VELOCITY)
                particles.append(
                    (
                        pos[0] + translate[0],
                        pos[1] + translate[1],
                        pos[2] + translate[2],
                        pos[3],
                        vel,
                    )
                )

        if include_elastic:
            return particles, copy.deepcopy(self.connections)
        else:
            return particles, []

    def add_particles(self, particles, connections):
        liquid_count = 0
        elast_count = 0
        boundary_count = 0

        for i in self.particles:
            if self.is_elastic(i):
                elast_count += 1
            elif self.is_boundary(i):
                boundary_count += 1
            elif self.is_liquid(i):
                liquid_count += 1

        for p in particles:
            self.particles[len(self.particles)] = {
                POSITION: (p[0], p[1], p[2], p[3]),
                VELOCITY: (p[4][0], p[4][1], p[4][2], p[4][3]),
            }
        for c in connections:
            # Connection IDs encode extra info in the fractional part (e.g. 1.3
            # means particle 1, sub-type .3).  NO_PARTICLE_ID is -1; any
            # non-negative integer part is a real particle reference that must
            # be offset.  Using int() handles 0.0 and 0.3 correctly; the old
            # `c[0] > 0` check incorrectly skipped particle 0.
            if int(c[0]) >= 0:
                c[0] = c[0] + elast_count

        self.connections += connections

    def __str__(self):
        pos_count = 0
        vel_count = 0

        liquid_count = 0
        elastic_count = 0
        boundary_count = 0
        type_counts = {}

        for i in self.particles:
            if POSITION in self.particles[i]:
                pos_count += 1
                type = self.particles[i][POSITION][3]
                if type not in type_counts:
                    type_counts[type] = 0
                type_counts[type] += 1

            if VELOCITY in self.particles[i]:
                vel_count += 1

            if self.is_liquid(i):
                liquid_count += 1
            elif self.is_elastic(i):
                elastic_count += 1
            elif self.is_boundary(i):
                boundary_count += 1
            else:
                raise ValueError(
                    f"Unknown particle type for particle {i}: {self.particles[i]}"
                )

        assert pos_count == vel_count, "Position and velocity counts do not match."
        assert pos_count == liquid_count + elastic_count + boundary_count, (
            "Particle counts do not match: "
            f"Positions: {pos_count}, Liquid: {liquid_count}, Elastic: {elastic_count}, Boundary: {boundary_count}"
        )

        info = (
            f"Configuration with {len(self.particles)} particles (liq: {liquid_count}, elast: {elastic_count}, bound: {boundary_count}), {len(self.connections)}(={len(self.connections) / 32}*32) connections, "
            f"{len(self.membranes)} membranes, and {len(self.particle_mem_index)} particle membrane indices."
        )
        info += f"\n    Simulation box: x: {self.simulation_box[0]}->{self.simulation_box[1]}, y: {self.simulation_box[2]}->{self.simulation_box[3]}, z: {self.simulation_box[4]}->{self.simulation_box[5]} "
        info += "\n    Type counts: %s" % dict(sorted(type_counts.items()))
        if self.physical_parameters:
            info += "\n    Physical parameters: %s" % self.physical_parameters

        return info

    def __repr__(self):
        return self.__str__()


def load_configuration_file(filename, verbose=False):
    configuration = Configuration()
    current_section = None

    print("Loading configuration from:", os.path.abspath(filename))
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Configuration file {filename} does not exist.")

    with open(filename, "r") as file:
        pos_count = 0
        vel_count = 0

        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Section headers are exactly "[name]".  Using startswith/endswith
            # avoids false-positives on data lines that happen to contain "[".
            if line.startswith("[") and line.endswith("]"):
                section_name = line[1:-1]
                if section_name == "end":
                    break
                current_section = section_name
                continue

            if current_section == SIMULATION_BOX:
                configuration.simulation_box.append(float(line))

            elif current_section == PHYSICAL_PARAMETERS:
                m = _PHYS_PARAM_RE.match(line)
                if m:
                    configuration.physical_parameters[m.group(1)] = float(m.group(2))

            elif current_section == POSITION:
                if verbose:
                    print("Checking position line: %s (%i)" % (line, pos_count))
                x, y, z, particle_type = line.split()
                configuration.particles[pos_count] = {
                    POSITION: (float(x), float(y), float(z), float(particle_type))
                }
                pos_count += 1

            elif current_section == VELOCITY:
                if verbose:
                    print("Checking velocity line: %s (%i)" % (line, vel_count))
                vx, vy, vz, particle_type = line.split()
                configuration.particles[vel_count][VELOCITY] = (
                    float(vx),
                    float(vy),
                    float(vz),
                    float(particle_type),
                )
                if verbose:
                    print(
                        f" Particle {vel_count} - pos: {configuration.particles[vel_count][POSITION]}; vel: {configuration.particles[vel_count][VELOCITY]}"
                    )
                vel_count += 1

            elif current_section == CONNECTION:
                a, b, c, d = line.split()
                configuration.connections.append(
                    [float(a), float(b), float(c), float(d)]
                )

            elif current_section == MEMBRANES:
                if verbose:
                    print("Checking membrane line: %s" % line)
                a, b, c = line.split()
                configuration.membranes.append((int(a), int(b), int(c)))

            elif current_section == PARTICLE_MEM_INDEX:
                if verbose:
                    print("Checking particle membrane index line: %s" % line)
                configuration.particle_mem_index.append(int(line))

    return configuration


def write_configuration_file(
    filename,
    configuration,
    include_velocity=True,
    include_liquid=True,
    include_elastics=True,
    include_boundary=True,
    verbose=False,
):
    with open(filename, "w") as file:
        if configuration.physical_parameters:
            file.write(f"[{PHYSICAL_PARAMETERS}]\n")
            for key, value in configuration.physical_parameters.items():
                file.write(f"{key}: {value}\n")

        file.write(f"[{SIMULATION_BOX}]\n")
        for value in configuration.simulation_box:
            file.write(f"{value}\n")

        file.write(f"[{POSITION}]\n")
        elastic_count = 0

        elast = ""
        elast_velocity = ""
        liquid = ""
        liquid_velocity = ""
        boundary = ""
        boundary_velocity = ""

        for i in configuration.particles:
            pos = configuration.particles[i][POSITION]
            p = f"{pos[0]}\t{pos[1]}\t{pos[2]}\t{pos[3]}\n"

            vel = configuration.particles[i][VELOCITY]
            v = f"{vel[0]}\t{vel[1]}\t{vel[2]}\t{vel[3]}\n"

            if configuration.is_elastic(i) and include_elastics:
                elast += p
                elast_velocity += v
            elif configuration.is_liquid(i) and include_liquid:
                liquid += p
                liquid_velocity += v
            elif configuration.is_boundary(i) and include_boundary:
                boundary += p
                boundary_velocity += v

            if configuration.is_elastic(i):
                elastic_count += 1

        file.write(elast)
        file.write(liquid)
        file.write(boundary)

        # Omit the [velocity] section entirely when include_velocity=False.
        # Writing the header with no entries would cause the C++ loader to read
        # the wrong number of velocity records.
        if include_velocity:
            file.write(f"[{VELOCITY}]\n")
            file.write(elast_velocity)
            file.write(liquid_velocity)
            file.write(boundary_velocity)

        file.write(f"[{CONNECTION}]\n")
        if include_elastics:
            # C++ allocates exactly numElasticP * MAX_NEIGHBOR_COUNT slots.
            assert (
                len(configuration.connections) == elastic_count * MAX_NEIGHBOR_COUNT
            ), (
                f"Expected {elastic_count * MAX_NEIGHBOR_COUNT} connections "
                f"({elastic_count} elastic × {MAX_NEIGHBOR_COUNT} slots), "
                f"got {len(configuration.connections)}."
            )
            for conn in configuration.connections:
                assert int(conn[0]) < elastic_count, (
                    f"Connection {conn} references particle {int(conn[0])} "
                    f"which meets or exceeds elastic count {elastic_count}."
                )
                file.write(f"{conn[0]}\t{conn[1]}\t{conn[2]}\t{conn[3]}\n")

        file.write(f"[{MEMBRANES}]\n")
        for mem in configuration.membranes:
            file.write(f"{mem[0]}\t{mem[1]}\t{mem[2]}\n")

        file.write(f"[{PARTICLE_MEM_INDEX}]\n")
        for index in configuration.particle_mem_index:
            file.write(f"{index}\n")

        file.write("[end]\n")


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 2 and sys.argv[1] == "-w2":
        config_file = "configuration/worm_alone_half_resolution"
        config_file = "configuration/worm_swim_half_resolution"

        out_file = "configuration/multiworm"
        conf = load_configuration_file(
            config_file,
        )
        conf2 = load_configuration_file(
            config_file,
        )

        print("-----")
        positions = [(25, 20, 5)]
        positions = [
            (25, 20, 5),
            (5, 30, 25),
            (5, 20, 65),
            (-15, 30, 65),
            (-5, 20, 5),
            (-25, 20, 65),
        ]
        positions = [
            (25, 20, 5),
            (5, 30, 25),
            (5, 20, 65),
        ]

        adds = []
        for p in positions:
            particles, connections = conf.get_particles(
                translate=p,
                include_liquid=False,
                include_elastic=True,
                include_boundary=False,
            )
            print(
                " Adding %s particles with %s connections to configuration."
                % (len(particles), len(connections))
            )
            conf2.add_particles(particles, connections)

        write_configuration_file(
            out_file,
            conf2,
            verbose=True,
            include_liquid=True,
            include_elastics=True,
            include_boundary=True,
        )
        print("Configuration loaded for test: %s" % conf)

        print("================================")
        conf3 = load_configuration_file(out_file)
        print("Configuration reloaded from output file: %s" % conf3)

    elif len(sys.argv) == 2 and sys.argv[1] == "-w":
        config_file_box = "configuration/demo1"
        config_file_worm = "configuration/worm_alone_half_resolution"
        out_file = "configuration/worm_box"

        conf_box = load_configuration_file(
            config_file_box,
        )
        conf_worm = load_configuration_file(
            config_file_worm,
        )

        print("-----")
        positions = []

        space = 20
        for x in range(-20, 60, space):
            for y in range(20, 120, space):
                for z in range(-20, 200, space):
                    print("Adding box at position: ", x, y, z)
                    positions.append((x, y, z))

        adds = []
        for p in positions:
            particles, connections = conf_box.get_particles(
                translate=p,
                include_liquid=True,
                include_elastic=False,
                include_boundary=False,
            )
            print(
                " Adding %s particles with %s connections to configuration."
                % (len(particles), len(connections))
            )
            conf_worm.add_particles(particles, connections)

        write_configuration_file(
            out_file,
            conf_worm,
            verbose=True,
            include_liquid=True,
            include_elastics=True,
            include_boundary=True,
        )
        print("Configuration written for worm test: %s" % conf_worm)

        print("================================")
        conf3 = load_configuration_file(out_file)
        print("Configuration reloaded from output file: %s" % conf3)

    elif len(sys.argv) == 2 and sys.argv[1] == "-g":
        config_file = "configuration/demo2"
        out_file = "configuration/gravity_test"
        conf = load_configuration_file(
            config_file,
        )
        write_configuration_file(
            out_file,
            conf,
            verbose=True,
            include_liquid=False,
            include_elastics=True,
            include_boundary=True,
        )
        print("Configuration loaded for gravity test: %s" % conf)

        print("================================")
        conf3 = load_configuration_file(out_file)
        print("Configuration reloaded from output file: %s" % conf3)

    else:
        if len(sys.argv) == 3:
            print("Usage: python ConfigBuilder.py <config_file> <output_file>")

            config_file = sys.argv[1]
            out_file = sys.argv[2]
            conf = load_configuration_file(
                config_file,
            )
            conf2 = load_configuration_file(
                config_file,
            )

            print("-----")
            poses = [(25, 20, 5), (5, 30, 25), (5, 20, 65)]

            adds = []
            for p in poses:
                particles, connections = conf.get_particles(
                    translate=p,
                    include_liquid=True,
                    include_elastic=True,
                    include_boundary=False,
                )
                print(
                    " Adding %s particles with %s connections to configuration."
                    % (len(particles), len(connections))
                )
                conf2.add_particles(particles, connections)

            """  """
            write_configuration_file(
                out_file,
                conf2,
                include_liquid=True,
                include_elastics=True,
                include_boundary=True,
            )
            conf3 = load_configuration_file(out_file)
            print("Configuration reloaded from output file: %s" % conf3)

        else:
            config_file = sys.argv[1]

            conf = load_configuration_file(
                config_file,
            )
            print("Configuration loaded: %s" % conf)
