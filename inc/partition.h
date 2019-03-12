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
            size_t end_cell_id;
            size_t start_ghost_cell_id;
            size_t end_ghost_cell_id;
            size_t size() const { return end - start; }
            size_t cell_count() const { return end_cell_id - start_cell_id; }
            size_t total_cell_count() const { return end_ghost_cell_id - start_ghost_cell_id; }
        };
    }
}

#endif //PROJECT_PARTITION_H
