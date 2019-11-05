import json


class Particle:
    def __init__(
            self,
            particle_id: int,
            cell_id: int,
    ):
        self.particle_id = particle_id
        self.cell_id = cell_id


def load_data(f_name):
    with open(f_name, 'r') as f:
        data = json.load(f)

        p_list = []
        for p in data:
            cell_id = p['cellId']
            p_id = p['particleId']
            p_list.append(Particle(cell_id=cell_id, particle_id=p_id))
        return p_list


def main() -> bool:
    s_p_list = load_data('s_thread')
    m_p_list = load_data('p_thread')
    bad = False
    for i in range(len(s_p_list)):
        if s_p_list[i].particle_id != m_p_list[i].particle_id:
            bad = True
            print(f"????? Bad  Example {i} ????")
            print(f"GOT p_id: {m_p_list[i].particle_id} cell_id {m_p_list[i].cell_id}. "
                  f"Expected p_id: {s_p_list[i].particle_id} cell_id:  {s_p_list[i].cell_id}")
            print("???????????????????????")
    if bad:
        return False

    return True


if __name__ == '__main__':
    print(main())
