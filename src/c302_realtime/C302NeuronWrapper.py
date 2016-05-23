# The MIT License (MIT)
#
# Copyright (c) 2011, 2013 OpenWorm.
# http://openworm.org
#
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the MIT License
# which accompanies this distribution, and is available at
# http://opensource.org/licenses/MIT
#
# Contributors:
#      OpenWorm - http://openworm.org/people.html
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
from __future__ import with_statement

import sys
import os.path

from NeuronWrapper import NrnSimulator
from neuron import gui  # TODO remove this than

from helper.myneuron import MyNeuron

__author__ = 'Sergey Khayrulin'

ca_levels = 'cai'

paramVec = [ca_levels]


class C302NrnSimulator(NrnSimulator):
    """
    NEURON wrapper specifically for the output of the c302 model. This is .py file that holds the calclium concentration
    of the muscle cells using parameters C1
    """

    def __init__(self, model_name, tstop=20.0, cell_names=None):
        if model_name != "":
            self.cell_names = cell_names
            if not (os.path.isfile(model_name)):
                raise AttributeError(
                    u"File: {0:s} doesn't exist please check the path to the file or name of file".format(model_name))
            # IN ORDER TO DO THIS YOU MUST SET YOUR PYTHONPATH TO YOUR nrn directory eg /usr/local/nrn/lib/python
            name = model_name[:-3]
            sys.path.insert(0, os.getcwd())
            module = __import__(name)

            h = module.h
            # print h
            # raise AttributeError()

            h.init()

            self.out_data = {}
            self.neurons_names = []
            self.neurons = {}
            self.sections = {}
            self.__find_all_neurons()

            if len(self.neurons_names) == 0:
                raise RuntimeError(u"In File: {0:s} with model no any neurons has been found. Please check the "
                                   u"the file".format(model_name))
            for name in self.neurons_names:  # TODO put check that we haven't added this neuron yet in dictionary neurons
                self.neurons[name] = MyNeuron(name, index=self.neurons_names.index(name))

            # Initialization of segments and data arrays
            for k, val in self.neurons.iteritems():
                val.init_sections(h, paramVec)
                # for sec in val.sections:
                #    self.sections[]
            self.__index_sub_segments()
            self.simulation_speed = 1



        else:
            raise ValueError("Name of file with Model shouldn't be empty")

    def __update_data(self):
        for k, val in self.neurons.iteritems():
            val.update_sec_data(paramVec)

    def one_step(self):
        """
        Make one step of NEURON simulation
        """
        from neuron import h
        if h.t < h.tstop:
            for i in xrange(self.simulation_speed):
                h.advance()
            self.__update_data()
            # gui.process_data()
        else:
            print 'Simulation is finished'


    def __find_all_neurons(self):
        """
        Search neurons names from hoc segment name
        """
        from neuron import h
        if(self.cell_names is not None):
            for cell in self.cell_names:
                var_name = "h.a_" + cell + "[0]"
                var = eval(var_name)
                self.neurons_names.append(var.hname())
        else:
            for h_sec in h.allsec():
                section_name = h_sec.name()
                self.neurons_names.append(section_name)

    def get_data(self):
        from neuron import h
        if h.t < h.tstop:
            data = []
            for n in self.neurons:
                data.append(self.neurons[n].get_averaged_data(ca_levels))
            return data
        else:
            return [0] * 95


    # def get_time(self):
    #     from neuron import h
    #     return h.t

    def __index_sub_segments(self):
        unique_indexes = []
        index = 0
        for k, v in self.neurons.iteritems():
            for sec in v.sections.values():
                for sub_sec in sec.sub_sections:
                    #if not(index in unique_indexes):
                    unique_indexes.append(index)
                    #else:
                    index += 1
                    sub_sec.index = index
        print index

    def get_max(self, cells):
        from neuron import h
        maxes = []
        while(h.t < h.tstop):
            self.one_step()
            l = max(self.get_data())
            maxes.append(l)
        h.t = 0
        return max(maxes)






    # def add_stim(self, n_name, sec_name=10):
    #     #TODO describe all in detail
    #     from neuron import h
    #     stim = h.IClamp(0.5, self.neurons[n_name].sections[sec_name].h_sec)
    #     stim.amp = 10.0
    #     stim.delay = 5.0
    #     stim.dur = 1.0

    # def add_synaps(self, sec_id1, sec_id2):
    #     pass

    # def finish(self):
    #     """
    #     Do nothing yet
    #     """
    #     pass
