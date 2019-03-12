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
            size_t start;
            size_t end;
            size_t start_cell_id;
            unsigned int end_cell_id;
            unsigned int start_ghost_cell_id;
            unsigned int end_ghost_cell_id;
            unsigned int size() const { return end - start; }
            unsigned int cell_count() const { return end_cell_id - start_cell_id; }
            unsigned int total_cell_count() const { return end_ghost_cell_id - start_ghost_cell_id; }
        };
    }
}

#endif //PROJECT_PARTITION_H
