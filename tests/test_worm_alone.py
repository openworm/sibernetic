from test_utils import plot_and_save_positions
import sys
import os


def test_cube1():

    sim_ref = "test_opencl_worm"
    types_to_plot = {1.1: [0, 100, 200, 300], 2.1: [0, 2]}
    plot_and_save_positions(
        sim_ref,
        sim_dir=os.path.abspath(f"simulations/{sim_ref}"),
        show_plot="-gui" in sys.argv,
        types_to_plot=types_to_plot,
    )


# main
if __name__ == "__main__":
    test_cube1()
