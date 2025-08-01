"""
Entry point for running the Sibernetic simulator alongside
c302 neuronal models. Provides command line options for configuring
duration, device selection and other parameters.
"""

try:
    from pyneuroml import pynml
except ImportError:  # pyneuroml may not be installed for simple tests
    pynml = None
import argparse
import re
import os
import sys
import time
import math

import pprint

pp = pprint.PrettyPrinter(indent=4)

script_version = "0.1.8"  # This will change at different rate to C++ code...

DEFAULTS = {
    "duration": 2.0,
    "dt": 0.005,
    "dtNrn": 0.05,
    "logstep": 1000,
    "reference": "Muscles",
    "c302params": "C1",
    "verbose": False,
    "device": "ALL",
    "configuration": "worm",
    "noc302": False,
    "outDir": "simulations",
    "datareader": "SpreadsheetDataReader",
    "test": False,
}


"""
    For recording the version of core Sibernetic in report file
"""


def get_sibernetic_version():
    # TODO make separate file for version info...
    version = None
    with open("src/main.cpp", "r") as sib_main_file:
        for line in sib_main_file:
            ws = line.split()
            if len(ws) >= 4 and ws[0] == "std::string" and ws[1] == "version":
                version = ws[3][1:-2]
                print_("Sibernetic v%s" % version)
    return version


def process_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "A script which can run Sibernetic, controlled by a "
            + "neuronal network in c302"
        )
    )

    parser.add_argument(
        "-test",
        action="store_true",
        default=DEFAULTS["test"],
        help="Use this flag to run some tests on the outputted files after the run, default: %s"
        % DEFAULTS["test"],
    )

    parser.add_argument(
        "-noc302",
        action="store_true",
        default=DEFAULTS["noc302"],
        help="Use this flag to run Sibernetic with the inbuilt (Python based) sine wave generator, using specified duration/configuration etc. and saving in simulations directory, default: %s"
        % DEFAULTS["noc302"],
    )

    parser.add_argument(
        "-duration",
        type=float,
        metavar="<duration>",
        default=DEFAULTS["duration"],
        help="Duration of simulation in ms, default: %sms" % DEFAULTS["duration"],
    )

    parser.add_argument(
        "-dt",
        type=float,
        metavar="<dt>",
        default=DEFAULTS["dt"],
        help="Time step for SIBERNETIC in ms, default: %sms" % DEFAULTS["dt"],
    )

    parser.add_argument(
        "-dtNrn",
        type=float,
        metavar="<dtNrn>",
        default=DEFAULTS["dtNrn"],
        help="Time step for NEURON in ms, default: %sms" % DEFAULTS["dtNrn"],
    )

    parser.add_argument(
        "-logstep",
        type=int,
        metavar="<logstep>",
        default=DEFAULTS["logstep"],
        help="Frequency of logging data into file, default: %s (note: Sibernetic's is 10)"
        % DEFAULTS["logstep"],
    )

    parser.add_argument(
        "-device",
        type=str,
        metavar="<device>",
        default=DEFAULTS["device"],
        help="Tell OpenCL to run on this device (could be cpu or gpu), default: %s"
        % DEFAULTS["device"],
    )

    parser.add_argument(
        "-configuration",
        type=str,
        metavar="<configuration>",
        default=DEFAULTS["configuration"],
        help="Run Sibernetic with this configuration (e.g. worm_no_water, worm_deep_water), default: %s"
        % DEFAULTS["configuration"],
    )

    parser.add_argument(
        "-reference",
        type=str,
        metavar="<reference>",
        default=DEFAULTS["reference"],
        help="Reference for c302 network subset (Muscles, Full, etc.), default: %s"
        % DEFAULTS["reference"],
    )

    parser.add_argument(
        "-c302params",
        type=str,
        metavar="<c302params>",
        default=DEFAULTS["c302params"],
        help="Parameter set from c302 (A, B, C, C1), default: %s"
        % DEFAULTS["c302params"],
    )

    parser.add_argument(
        "-datareader",
        type=str,
        metavar="<datareader>",
        default=DEFAULTS["datareader"],
        help="DataReader, default: %s" % DEFAULTS["datareader"],
    )

    parser.add_argument(
        "-outDir",
        type=str,
        metavar="<outDir>",
        default=DEFAULTS["outDir"],
        help="Output directory, default: %s" % DEFAULTS["outDir"],
    )

    return parser.parse_args()


def print_(msg):
    pre = "Sib_c302  >>>"
    print("%s %s" % (pre, msg.replace("\n", "\n" + pre + " ")))


def main(args=None):
    if args is None:
        args = process_args()
    run(a=args)


def build_namespace(a=None, **kwargs):
    if a is None:
        a = argparse.Namespace()

    # Add arguments passed in by keyword.
    for key, value in kwargs.items():
        setattr(a, key, value)

    # Add defaults for arguments not provided.
    for key, value in DEFAULTS.items():
        if not hasattr(a, key):
            setattr(a, key, value)

    # Change all values to under_score from camelCase.
    name_changes = []
    for key, value in a.__dict__.items():
        new_key = convert_case(key)
        if new_key != key:
            name_changes.append((key, new_key, value))
    for old_key, new_key, value in name_changes:
        setattr(a, new_key, value)
        delattr(a, old_key)

    return a


def convert_case(name):
    """Converts from camelCase to under_score"""
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


def announce(message):
    print_(
        "\n************************************************************************\n*"
    )
    print_("*  %s" % message.replace("\n", "\n*  "))
    print_(
        "*\n************************************************************************"
    )


def check_file(file_name, dir, line_num=None):
    file_name = os.path.join(dir, file_name)
    print_(
        "Checking %s exists %s"
        % (file_name, "(and has %s lines)" % line_num if line_num else "")
    )

    if not os.path.isfile(file_name):
        print_("Error: File %s missing!!" % file_name)
        return False
    if not os.path.getsize(file_name) > 0:
        print_("Error: File %s is empty!!" % file_name)
        return False
    if line_num:
        count = 0
        for line in open(file_name):
            count += 1
        if count != line_num:
            print_("Error: Expecting %s lines, but had %s!!" % (line_num, count))
            return False

    return True


def dynamic_import(abs_module_path, class_name):
    from importlib import import_module

    module_object = import_module(abs_module_path)

    target_class = getattr(module_object, class_name)

    return target_class


def run(a=None, **kwargs):
    a = build_namespace(a, **kwargs)

    if not a.noc302:
        try:
            import neuroml  # noqa: F401
            import pyneuroml  # noqa: F401
            import xlrd  # noqa: F401
        except Exception as e:
            print_(
                "Cannot import one of the required packages. Please install!\n"
                "Exception: %s\n" % e
            )

    if not a.noc302:
        try:
            if "C302_HOME" in os.environ:
                os.environ["C302_HOME"]
                sys.path.append(os.environ["C302_HOME"])
                print_("Python path now: %s" % sys.path)
            import c302
            import c302.c302_utils as c302_utils
        except Exception as e:
            print_(
                "Cannot import c302!\n"
                "Exception: %s\n"
                % e
                + "Please see installation instructions at https://github.com/openworm/c302!\n"
            )

            exit(-1)

    gen_start = time.time()

    # if not os.path.isdir('simulations'):
    #    os.mkdir('simulations')
    if not os.path.isdir(a.out_dir):
        os.mkdir(a.out_dir)

    current_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.gmtime())

    if not a.noc302:
        ref = a.reference
        sim_ref = "%s_%s_%s" % (a.c302params, ref, current_time)

    else:
        sim_ref = "Sibernetic_%s" % (current_time)

    # sim_dir = "simulations/%s" % (sim_ref)
    sim_dir = os.path.join(a.out_dir, sim_ref)

    os.makedirs(sim_dir, exist_ok=True)

    run_dir = "."
    if "SIBERNETIC_HOME" in os.environ:
        run_dir = os.environ["SIBERNETIC_HOME"]

    if not a.noc302:
        id = "%s_%s" % (a.c302params, ref)

        setup = dynamic_import("c302.c302_%s" % ref, "setup")

        setup(
            a.c302params,
            generate=True,
            duration=a.duration,
            dt=a.dt,
            target_directory=sim_dir,
            data_reader=a.datareader,
        )

        lems_file0 = os.path.join(sim_dir, "LEMS_c302_%s.xml" % id)
        lems_file = os.path.join(sim_dir, "LEMS_c302.xml")
        print_("Renaming %s -> %s" % (lems_file0, lems_file))
        os.rename(lems_file0, lems_file)

        announce("Generating NEURON files from: %s..." % lems_file)

        pynml.run_lems_with_jneuroml_neuron(
            lems_file,
            only_generate_scripts=True,
            nogui=True,
            load_saved_data=False,
            verbose=True,
            realtime_output=True,
        )

        """with open(os.path.join(sim_dir, 'LEMS_c302_nrn.py'), 'r') as main_nrn_py:
            updated =''
            for line in main_nrn_py:
                line = line.replace('GenericCell.hoc','%s/GenericCell.hoc' % sim_dir)
                line = line.replace('GenericNeuronCell.hoc','%s/GenericNeuronCell.hoc' % sim_dir)
                line = line.replace('GenericMuscleCell.hoc','%s/GenericMuscleCell.hoc' % sim_dir)
                line = line.replace("open('time.dat","open('%s/time.dat" % sim_dir)
                line = line.replace("open('c302_","open('%s/c302_" % sim_dir)
                updated += line"""
        main_nrn_py = open("%s/LEMS_c302_nrn.py" % (sim_dir), "r")
        updated = ""
        for line in main_nrn_py:
            line = line.replace("GenericCell.hoc", "%s/GenericCell.hoc" % sim_dir)
            line = line.replace(
                "GenericNeuronCell.hoc", "%s/GenericNeuronCell.hoc" % sim_dir
            )
            line = line.replace(
                "GenericMuscleCell.hoc", "%s/GenericMuscleCell.hoc" % sim_dir
            )
            line = line.replace("open('time.dat", "open('%s/time.dat" % sim_dir)
            line = line.replace("open('c302_", "open('%s/c302_" % sim_dir)
            updated += line
        main_nrn_py.close()

        """with open(os.path.join(sim_dir, 'LEMS_c302_nrn.py'), 'w') as main_nrn_py:
            main_nrn_py.write(updated)"""
        main_nrn_py = open("%s/LEMS_c302_nrn.py" % (sim_dir), "w")
        main_nrn_py.write(updated)
        main_nrn_py.close()

        command = "nrnivmodl %s" % sim_dir
        # command = 'nrnivmodl .'
        announce("Compiling NMODL files for NEURON...")
        try:
            pynml.execute_command_in_dir_with_realtime_output(
                command, run_dir, prefix="nrnivmodl >> "
            )
        except KeyboardInterrupt:
            print_("\nCaught CTRL+C\n")
            sys.exit()

    command = (
        "./Release/Sibernetic %s -f %s -no_g -l_to lpath=%s timelimit=%s timestep=%s logstep=%s device=%s"
        % (
            "" if a.noc302 else "-c302",
            a.configuration,
            sim_dir,
            a.duration / 1000.0,
            a.dt / 1000,
            a.logstep,
            a.device,
        )
    )

    env = {
        "DISPLAY": os.environ.get("DISPLAY") if os.environ.get("DISPLAY") else "",
        "XAUTHORITY": os.environ.get("XAUTHORITY")
        if os.environ.get("XAUTHORITY")
        else "",
        "PYTHONPATH": ".:%s:%s"
        % (os.environ.get("PYTHONPATH", "."), os.path.abspath(sim_dir)),
        "NEURON_MODULE_OPTIONS": "-nogui",
        "PATH": os.environ.get("PATH", "/usr/bin"),
    }

    sim_start = time.time()

    reportj = {}

    SUCCESS = "Completed successfully"
    completion_status = "Not completed"

    announce(
        "Executing main Sibernetic simulation of %sms using: \n\n    %s \n\n  in %s with %s"
        % (a.duration, command, run_dir, env)
    )
    try:
        pynml_success = pynml.execute_command_in_dir_with_realtime_output(
            command, run_dir, prefix="Sibernetic >> ", env=env, verbose=True
        )
        completion_status = (
            SUCCESS if pynml_success else "Failure running command: %s" % command
        )
    except KeyboardInterrupt:
        print_("\nCaught CTRL+C. Continue...\n")
        completion_status = "Caught CTRL+C"
    except Exception as e:
        print_("\nError in running Sibernetic: %s...\n" % e)
        completion_status = "Error during simulation"

    reportj["completion_status"] = completion_status

    successful = completion_status == SUCCESS

    sim_end = time.time()

    reportj["duration"] = "%s ms" % a.duration
    reportj["dt"] = "%s ms" % a.dt
    reportj["dtNrn"] = "%s ms" % a.dt_nrn
    reportj["logstep"] = a.logstep
    reportj["sim_ref"] = sim_ref
    reportj["noc302"] = a.noc302
    reportj["python_args"] = "python " + " ".join(sys.argv)
    reportj["sibernetic_version"] = get_sibernetic_version()
    reportj["sibernetic_c302_version"] = script_version
    import platform

    reportj["python_version"] = platform.python_version()
    reportj["os_version"] = ", ".join(platform.uname())

    if not a.noc302:
        reportj["reference"] = a.reference
        reportj["c302params"] = a.c302params
        reportj["c302_version"] = c302.__version__

        for m in ["pyneuroml", "neuroml", "matplotlib", "numpy", "cect"]:
            if m == "neuroml":
                m_ = "libneuroml"
            else:
                m_ = m
            try:
                exec("import %s" % m)
                installed_ver = "%s" % eval("%s.__version__" % m)
                reportj["%s_version" % m_] = installed_ver
            except Exception:
                reportj["%s_version" % m_] = "-- Not found! --"

        from neuron import h

        reportj["neuron_version"] = h.nrnversion()

        try:
            import PyOpenWorm as P

            reportj["pyopenworm_version"] = P.__version__
        except Exception:
            pass  # Not currently required for sims

    reportj["generation_time"] = "%s s" % (sim_start - gen_start)
    reportj["run_time"] = "%s s" % (sim_end - sim_start)
    reportj["device"] = "%s" % (a.device)
    reportj["configuration"] = "%s" % (a.configuration)
    reportj["run_command"] = command
    reportj["outDir"] = a.out_dir
    reportj["env"] = str(env)

    with open(os.path.join(sim_dir, "report.json"), "w") as report_file:
        import json

        s = json.dumps(reportj, indent=4, sort_keys=True)
        report_file.write(s)

    if not a.noc302 and successful:
        announce(
            "Generating images for neuronal activity (via %s in %s)..."
            % (lems_file, sim_dir)
        )

        results = pynml.reload_saved_data(
            lems_file,
            base_dir=sim_dir,
            plot=False,
            show_plot_already=False,
            simulator=None,
            verbose=True,
        )

        c302_utils.plot_c302_results(
            results,
            config=a.reference,
            parameter_set=a.c302params,
            directory=sim_dir,
            save=True,
            show_plot_already=False,
        )

    ####   WCON plotting better...
    # pos_file_name = os.path.abspath('%s/position_buffer.txt'%sim_dir)
    # announce("Plotting positions of worm body particles in %s..."%pos_file_name)
    # from plot_positions import plot_positions
    # if not os.path.isfile(pos_file_name):
    #    time.sleep(2)
    # plot_positions(pos_file_name,rate_to_plot = int(a.duration/5), show_plot=False)

    if successful:
        from plot_positions import plot_muscle_activity

        musc_act_file = os.path.join(sim_dir, "muscles_activity_buffer.txt")
        plot_muscle_activity(musc_act_file, a.dt, a.logstep, show_plot=False)

        from wcon.generate_wcon import generate_wcon

        num_steps = int(a.duration / a.dt)
        num_steps_logged = num_steps / a.logstep
        rate = max(1, int(num_steps_logged / 20.0))
        # print("%s, %s"%(num_steps, rate))
        generate_wcon(
            os.path.join(sim_dir, "worm_motion_log.txt"),
            os.path.join(sim_dir, "worm_motion_log.wcon"),
            rate_to_plot=rate,
            plot=False,
            save_figure1_to=os.path.join(sim_dir, "worm_motion_1.png"),
            save_figure2_to=os.path.join(sim_dir, "worm_motion_2.png"),
            save_figure3_to=os.path.join(sim_dir, "worm_motion_3.png"),
        )

    run_dur_sec = sim_end - sim_start
    announce(
        "Finished run in %s sec (%s hours)!\n\n" % (run_dur_sec, run_dur_sec / 3600.0)
        + (
            "Simulation output saved in: %s\n\n" % (sim_dir)
            if successful
            else "SIMULATION UNSUCCESSFUL!\n\n"
        )
        + "Report of simulation at: %s/report.json\n\n" % (sim_dir)
        + (
            "Replay recorded simulation with: ./Release/Sibernetic -l_from lpath=%s\n"
            % (sim_dir)
            if successful
            else ""
        )
    )

    if a.test:
        announce("Running tests...")
        passed = True

        total_steps = int(a.duration / a.dt)
        sib_recorded_steps = int(math.ceil(float(total_steps) / a.logstep))

        print_(
            "Expected Sibernetic to have recorded %s steps of %s total in %s ms..."
            % (sib_recorded_steps, total_steps, a.duration)
        )

        passed = passed and check_file(
            "worm_motion_log.txt", sim_dir, line_num=sib_recorded_steps
        )
        passed = passed and check_file("worm_motion_log.wcon", sim_dir)
        passed = passed and check_file("worm_motion_1.png", sim_dir)
        passed = passed and check_file("worm_motion_2.png", sim_dir)
        passed = passed and check_file("worm_motion_3.png", sim_dir)

        passed = passed and check_file("connection_buffer.txt", sim_dir)
        passed = passed and check_file("membranes_buffer.txt", sim_dir)
        passed = passed and check_file(
            "muscles_activity_buffer.txt", sim_dir, line_num=sib_recorded_steps
        )
        passed = passed and check_file("position_buffer.txt", sim_dir)
        passed = passed and check_file("report.json", sim_dir)

        if not passed:
            announce("Failed tests!!")
            exit(-1)
        else:
            announce("Passed all tests!!")


if __name__ == "__main__":
    main()
