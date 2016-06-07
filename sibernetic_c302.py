
from pyneuroml import pynml
import argparse

DEFAULTS = {'duration': 100,
            'verbose': False} 
            
def process_args():
    """ 
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
                description=("A script which can be run to generate a LEMS "
                             "file to analyse the behaviour of channels in "
                             "NeuroML 2"))

    parser.add_argument('-duration', 
                        type=float,
                        metavar='<duration>',
                        default=DEFAULTS['duration'],
                        help="Duration of simulation in ms, default: %sms"%DEFAULTS['duration'])
                        

command = './Release/Sibernetic -f worm timelimit=0.001'

pynml.execute_command_in_dir(command, '.', prefix="Sibernetic: ")




