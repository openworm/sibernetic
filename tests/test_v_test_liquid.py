from test_utils import plot_and_save_positions
import sys
import os


def test_v_test_liquid():

    sim_ref = "test_opencl_v_test_liquid"
    types_to_plot = {1.1: [1, 200, 72100, 72000]}
    plot_and_save_positions(
        sim_ref,
        sim_dir=os.path.abspath(f"../simulations/{sim_ref}"),
        show_plot="-gui" in sys.argv,
        types_to_plot=types_to_plot,
        index_to_save=1,
    )


# main
if __name__ == "__main__":
    test_v_test_liquid()
