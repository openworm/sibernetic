import json
import math

from typing import List

NEIGBOUR_COUNT = 32
H = 3.34
H_2 = H ** 2
SIM_SCALE = 7.4e-06
H_SCALED = H * SIM_SCALE 

EPS = 0.000001

class Float3:
    def __init__(
        self,
        x,
        y,
        z
    ) -> None:
        self.x = x
        self.y = y
        self.z = z
    
    def dot(self):
        return self.x * self.x + self.y * self.y + self.z * self.z

    @staticmethod
    def dist(p1, p2):
        tmp = Float3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
        return math.sqrt(tmp.dot())
    
    @staticmethod
    def dist_2(p1, p2):
        tmp = Float3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
        return tmp.dot()

class Neighbour:
    def __init__(
        self,
        n_particle_id: float,
        distance: float
    ) -> None:
        self.p_id = n_particle_id
        self.dist = distance
    
    def __eq__(self, n) -> bool:
        return self.p_id == n.p_id and abs(self.dist - n.dist) <= EPS

class Particle:
    def __init__(
        self, 
        json_particle
    ) -> None:
        self.pos = Float3(**json_particle['particle']['pos'])
        self.vel = Float3(**json_particle['particle']['vel'])
        self.accel = Float3(**json_particle['particle']['accel'])
        self.density = json_particle['particle']['density']
        self.pressure = json_particle['particle']['pressure']
        self.viscosity = json_particle['particle']['viscosity']
        self.mass = json_particle['particle']['mass']
        self.cell_id = json_particle['particle']['cell_id']
        # if self.density > 1000.0:
        #     print(self.density)
        #     print(self.pressure)
        # if self.mass > 1.9999999920083944E-12:
        #     print(self.mass)
        self.id = json_particle['particle']['particle_id']
        self.stored_n_map = {}
        for json_n in json_particle['n_list']:
            p_id = int(json_n['n_particle_id'])
            if p_id in self.stored_n_map:
                self.stored_n_map[p_id].append(Neighbour(**json_n))
            
            else:
                self.stored_n_map[p_id] = [Neighbour(**json_n)]

    def calc_n_map(self, particle_list: List['Particle']) -> None:
        tmp_n_list = {}
        p_id = 0
        for p in particle_list:
            if p != self:
                d = Float3.dist_2(self.pos, p.pos)
                if d <= H_2:
                    tmp_n_list[p_id] = Neighbour(n_particle_id=p_id, distance=math.sqrt(d) * SIM_SCALE)
            p_id += 1
        if len(tmp_n_list) > NEIGBOUR_COUNT:
            tmp_dct = [(k, v) for k, v in tmp_n_list.items()]
            tmp_dct.sort(key=lambda p: p[1].dist)
            tmp_n_list = dict(tmp_dct[:NEIGBOUR_COUNT])
        elif len(tmp_n_list) < NEIGBOUR_COUNT:
            tmp_n_list[-1] = [
                Neighbour(-1.0, -1.0)
                for _ in range(NEIGBOUR_COUNT - len(tmp_n_list))
            ]
        #tmp_n_list.sort(key=lambda p: p.p_id, reverse=True)
        self.n_map = tmp_n_list

    def get_p(self, real_p_id, particle_list):
        for p in particle_list:
            if p.id == real_p_id:
                return p

    def compare_neighbour_list(self, particle_list: List['Particle']) -> bool:
        bad = True
        diff = set(self.n_map.keys()).difference(self.stored_n_map.keys())

        if diff:
            print("DIff is", diff)
            print("K1", self.n_map.keys())
            print("K2", self.stored_n_map.keys())
            return False
        return bad

def test_nmap(particle_list: List[Particle]) -> bool:
    i = 0
    result = True
    for p in particle_list:
        if p.compare_neighbour_list(particle_list) == False:
            print("indx ", i)
            result = False
        i += 1

    return result

def test_nmap_particle(p, particle_list):
    if p.compare_neighbour_list(particle_list) == False:
        return False

    return True

def main() -> None:
    ns_map = None
    with open("../../debug", 'r') as f:
        ns_map = json.load(f)
    if ns_map is not None:
        particle_list = [Particle(json_p) for json_p in ns_map]
        particle_to_check = len(particle_list)
        _all = True
        if _all:
            for p in particle_list:
                p.calc_n_map(particle_list)
            print(test_nmap(particle_list))
        else:
            p = particle_list[34]
            p.calc_n_map(particle_list)
            print(test_nmap_particle(p, particle_list))
        
        print("Checked ", particle_to_check)
    else:
        print(False)

    

if __name__ == '__main__':
    main()