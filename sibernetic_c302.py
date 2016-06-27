from pyneuroml import pynml
import argparse
import re
import os
import sys
import time

DEFAULTS = {'duration': 1.0,
            'dt': 0.005,
            'dtNrn': 0.05,
            'reference': 'Muscles',
            'c302params': 'C1',
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
                        help="Reference for network subset (Muscles, Full, etc.), default: %s"%DEFAULTS['reference'])
                        
    parser.add_argument('-c302params', 
                        type=str,
                        metavar='<c302params>',
                        default=DEFAULTS['c302params'],
                        help="Parameter set from c302 (A, B, C, C1), default: %s"%DEFAULTS['c302params'])
                        
    return parser.parse_args()
                        



def print_(msg):
    pre = "Sib_c302  >>>"
    print('%s %s'%(pre,msg.replace('\n','\n'+pre+' ')))


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
    
    try:
        if os.environ.has_key('C302_HOME'):
            os.environ['C302_HOME']
            sys.path.append(os.environ['C302_HOME'])
            print_('Python path now: %s'%sys.path)
        import c302
    except:
        print_("Cannot import c302!\n"
             +"Please set environment variable C302_HOME to point to the directory: CElegansNeuroML/CElegans/pythonScripts/c302!\n")
             
        exit()
        
    a = build_namespace(a,**kwargs)
    
    
    ref = a.reference
    
    if not os.path.isdir('simulations'):
        os.mkdir('simulations')
    
    sim_dir = "simulations/%s_%s"%(ref, time.ctime().replace(' ','_' ).replace(':','.' ))
    os.mkdir(sim_dir)
    
    #exec('from %s import ParameterisedModel'%a.c302params)
    #params = ParameterisedModel()
    
    id = '%s_%s'%(a.c302params,ref)
    
    
    exec('from c302_%s import setup'%ref)
    
    setup(a.c302params, 
          generate=True,
          duration = a.duration,
          dt = a.dt,
          target_directory=sim_dir)
    
             
    lems_file0 = '%s/LEMS_c302_%s.xml'%(sim_dir,id)
    lems_file = '%s/LEMS_c302.xml'%(sim_dir)
    print_("Renaming %s -> %s"%(lems_file0,lems_file))
    os.rename(lems_file0,lems_file)
    
    print_("Generating NEURON files from: %s"%lems_file)
    
    pynml.run_lems_with_jneuroml_neuron(lems_file,
                                        only_generate_scripts=True,
                                        nogui=True, 
                                        load_saved_data=False, 
                                        verbose=True)
                                        
    main_nrn_py = open('%s/LEMS_c302_nrn.py'%(sim_dir),'r')
    updated =''
    for line in main_nrn_py:
        updated += '%s'%line.replace('GenericCell.hoc','%s/GenericCell.hoc'%sim_dir)
    main_nrn_py.close() 
    
    main_nrn_py = open('%s/LEMS_c302_nrn.py'%(sim_dir),'w')
    main_nrn_py.write(updated)
    main_nrn_py.close() 
    
    run_dir = '.'
    command = 'nrnivmodl %s'%sim_dir

    print_("Executing: %s in %s"%(command, ''))
    pynml.execute_command_in_dir(command, run_dir, prefix="nrnivmodl: ")

    command = './Release/Sibernetic -f worm -no_g -l_to lpath=%s timelimit=%s timestep=%s'%(sim_dir,a.duration/1000.0,a.dt/1000)
    env={"PYTHONPATH":"./src:./%s"%sim_dir}
    print_("Executing: %s in %s with %s"%(command, run_dir, env))
    #pynml.execute_command_in_dir('env', run_dir, prefix="Sibernetic: ",env=env,verbose=True)
    pynml.execute_command_in_dir(command, run_dir, prefix="Sibernetic: ",env=env,verbose=True)
    
    
    print_("Finished!|\nRerun simulation with: ./Release/Sibernetic -l_from lpath=%s\n"%(sim_dir))



if __name__ == '__main__':
    main()


