from pyneuroml import pynml
import argparse
import re
import os
import time

import c302

DEFAULTS = {'duration': 100.0,
            'dt': 0.005,
            'dtNrn': 0.05,
            'reference': 'Muscles',
            'c302params': 'parameters_C1',
            'verbose': False} 
            
            
def process_args():
    """ 
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
                description=("A script which can run Sibernetic, controlled by a "+
                "neuronal network in c302"))

    parser.add_argument('-duration', 
                        type=float,
                        metavar='<duration>',
                        default=DEFAULTS['duration'],
                        help="Duration of simulation in ms, default: %sms"%DEFAULTS['duration'])
                        
    parser.add_argument('-dt', 
                        type=float,
                        metavar='<dt>',
                        default=DEFAULTS['dt'],
                        help="Time step for SIBERNETIC in ms, default: %sms"%DEFAULTS['dt'])
                        
    parser.add_argument('-dtNrn', 
                        type=float,
                        metavar='<dtNrn>',
                        default=DEFAULTS['dtNrn'],
                        help="Time step for NEURON in ms, default: %sms"%DEFAULTS['dtNrn'])
                        
    parser.add_argument('-reference', 
                        type=str,
                        metavar='<reference>',
                        default=DEFAULTS['reference'],
                        help="String to use in simulation results directory, default: %s"%DEFAULTS['reference'])
                        
    parser.add_argument('-c302params', 
                        type=str,
                        metavar='<c302params>',
                        default=DEFAULTS['c302params'],
                        help="Parameter set from c302, default: %s"%DEFAULTS['c302params'])
                        
    return parser.parse_args()
                        


def print_(msg):
    pre = "Sib_c302 >> "
    print('%s %s'%(pre,msg.replace('\n','\n'+pre)))


def main(args=None):
    if args is None:
        args = process_args()
    run(a=args)
    

def build_namespace(a=None,**kwargs):
    if a is None:
        a = argparse.Namespace()
    
    # Add arguments passed in by keyword.  
    for key,value in kwargs.items():
        setattr(a,key,value)

    # Add defaults for arguments not provided.  
    for key,value in DEFAULTS.items():
        if not hasattr(a,key):
            setattr(a,key,value)

    # Change all values to under_score from camelCase.  
    for key,value in a.__dict__.items():
        new_key = convert_case(key)
        if new_key != key:
            setattr(a,new_key,value)
            delattr(a,key)

    return a


def convert_case(name):
    """Converts from camelCase to under_score"""
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def run(a=None,**kwargs): 
    a = build_namespace(a,**kwargs)
    
    
    ref = a.reference
    
    run_dir = "simulations/%s_%s"%(ref, time.ctime().replace(' ','_' ).replace(':','.' ))
    os.mkdir(run_dir)
    run_dir = "simulations/%s_%s"%(ref, 0)
    
    exec('from %s import ParameterisedModel'%a.c302params)
    params = ParameterisedModel()
    
    id = '%s_%s'%(ref,a.c302params.split('_')[1])
    
    c302.generate(id,
             params,
             cells = None,
             cells_to_plot = None,
             cells_to_stimulate = None,
             include_muscles=True,
             conn_number_override = None,
             conn_number_scaling = None,
             duration = a.duration,
             dt = a.dt,
             seed = 1234,
             validate=True, 
             test=False,
             verbose=True,
             target_directory=run_dir)
             
    lems_file0 = '%s/LEMS_%s.xml'%(run_dir,id)
    lems_file = '%s/LEMS_c302.xml'%(run_dir)
    os.rename(lems_file0,lems_file)
    
    print_("Generating NEURON files from: %s"%lems_file)
    
    pynml.run_lems_with_jneuroml_neuron(lems_file,
                                        only_generate_scripts=True,
                                        nogui=True, 
                                        load_saved_data=False, 
                                        verbose=True)
    

    command = 'nrnivmodl'

    print_("Executing: %s in %s"%(command, run_dir))
    pynml.execute_command_in_dir(command, run_dir, prefix="nrnivmodl: ")

    command = '../../Release/Sibernetic -f worm timelimit=%s timestep=%s'%(a.duration/1000.0,a.dt/1000)

    print_("Executing: %s in %s"%(command, run_dir))
    #pynml.execute_command_in_dir(command, run_dir, prefix="Sibernetic: ",verbose=True)



if __name__ == '__main__':
    main()


