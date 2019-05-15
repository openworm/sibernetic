import json
import math

from typing import List

NEIGBOUR_COUNT = 32
H = 3.34
H_2 = H ** 2
SIM_SCALE = 0.06508344313069617
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

        self.id = json_particle['particle_id']
        self.stored_n_map = [
            Neighbour(**json_n)
            for json_n in json_particle['n_list']
        ]
        self.stored_n_map.sort(key=lambda p: p.p_id)

    def calc_n_map(self, particle_list: List['Particle']) -> None:
        tmp_n_list = []
        for p in particle_list:
            if p != self:
                d = Float3.dist_2(self.pos, p.pos)
                if d <= H_2:
                    tmp_n_list.append(Neighbour(n_particle_id=p.id, distance=math.sqrt(d) * SIM_SCALE))
        
        if len(tmp_n_list) > NEIGBOUR_COUNT:
            tmp_n_list = tmp_n_list[:NEIGBOUR_COUNT]
        else:
            tmp_n_list += [
                Neighbour(-1.0, -1.0)
                for _ in range(NEIGBOUR_COUNT - len(tmp_n_list))
            ]
        tmp_n_list.sort(key=lambda p: p.p_id)
        self.n_map = tmp_n_list

    def compare_neighbour_list(self) -> bool:
        for i, n in enumerate(self.stored_n_map):
            if n != self.n_map[i]:
                print(f"particle id - {self.id}, incorrect neigh {i}" )
                #print(f"{}")
                return False

        return True

def test_nmap(particle_list: List[Particle]) -> bool:
    for p in particle_list:
        if p.compare_neighbour_list() == False:
            return False

    return True

def main() -> None:
    ns_map = None
    with open("../../debug", 'r') as f:
        ns_map = json.load(f)
    if ns_map is not None:
        particle_list = [Particle(json_p) for json_p in ns_map]
        particle_to_check = len(particle_list)#500
        for p in particle_list[:particle_to_check]:
            p.calc_n_map(particle_list)
        
        print(test_nmap(particle_list[:particle_to_check]))
        print("Checked ", particle_to_check)
    else:
        print(False)

    

if __name__ == '__main__':
    main()