#include <iostream>
#include <memory>
#include <vector>
#include "ocl_helper.h"



int main(int argc, char **argv){
    p_q dev_q = get_dev_queue();
    std::vector<std::shared_ptr<i_solver>> _solvers;
    while (!dev_q.empty()) {
        try {
            std::shared_ptr<ocl_solver<T>> solver(
                    new ocl_solver<T>(model, dev_q.top(), device_index));
            _solvers.push_back(solver);
            std::cout << "************* DEVICE *************" << std::endl;
            dev_q.top()->show_info();
            std::cout << "**********************************" << std::endl;
        } catch (ocl_error &ex) {
            std::cout << ex.what() << std::endl;
        }
        dev_q.pop();
    }
    return 0;
}