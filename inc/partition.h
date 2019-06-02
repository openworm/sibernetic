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
            int start_cell_id;
            int end_cell_id;
            int start_ghost_cell_id;
            int end_ghost_cell_id;
            int size() const {
            	if(start == 0)
            	    return end - start;
	            return end - start + 1;
            }
            int total_size() const {
	            if(ghost_start == 0)
		            return (ghost_end - ghost_start);
            	return (ghost_end - ghost_start + 1);
            }
            uint cell_count() const { return end_cell_id - start_cell_id + 1; }
            int total_cell_count() const {
	            if(start_ghost_cell_id == 0)
	                return end_ghost_cell_id - start_ghost_cell_id;
	            return end_ghost_cell_id - start_ghost_cell_id;// + 1;
            }
            int offset() {
                return static_cast<int>((start > ghost_start) ? (start - ghost_start) : 0);
            }
            int limit() {
                return size() + offset();
            }
        };
    }
}

#endif //PROJECT_PARTITION_H
