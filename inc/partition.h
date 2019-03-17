//
// Created by sergey on 13.02.19.
//

#ifndef PROJECT_PARTITION_H
#define PROJECT_PARTITION_H

#include <iostream>

namespace sibernetic{
    namespace model{
        struct partition {
            /**each device has its own partition
             * in which we define where starts
             * and end particles for this device.
             */
            typedef unsigned int uint;
            size_t start;
            size_t end;
            size_t ghost_start;
            size_t ghost_end;
            uint start_cell_id;
            uint end_cell_id;
            uint start_ghost_cell_id;
            uint end_ghost_cell_id;
            uint size() const { return static_cast<uint>(end - start); }
            uint total_size() const { return static_cast<uint>(ghost_end - ghost_start ); }
            uint cell_count() const { return end_cell_id - start_cell_id + 1; }
            uint total_cell_count() const { return end_ghost_cell_id - start_ghost_cell_id + 1; }
        };
    }
}

#endif //PROJECT_PARTITION_H
