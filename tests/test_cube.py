from test_utils import plot_and_save_positions
import sys
import os


def test_cube():

    sim_ref = "test_cube"
    types_to_plot = {1.1: [0, 60, 124], 2.1: [0, 60, 124]}
    plot_and_save_positions(
        sim_ref,
        sim_dir=os.path.abspath(f"simulations/{sim_ref}"),
        show_plot="-gui" in sys.argv,
        types_to_plot=types_to_plot,
    )


# main
if __name__ == "__main__":
    test_cube()
