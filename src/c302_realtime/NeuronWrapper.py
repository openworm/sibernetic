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


from neuron import gui #TODO remove this than
from helper.myneuron import MyNeuron

__author__ = 'Sergey Khayrulin'

v = 'v'
v_pre = 'v_pre'
v_post = 'v_post'
i_syn = 'i_syn'
t = 't'
paramVec = [v]



class NrnSimulator:
    """
    NEURON wrapper init and run
    model described on hoc into NEURON simulator
    """

    def __init__(self, model_name, tstop=20.0):
        if model_name != "":
            if not (os.path.isfile(model_name)):
                raise AttributeError(
                    u"File: {0:s} doesn't exist please check the path to the file or name of file".format(model_name))
            from neuron import h
            h.finitialize()
            h.load_file("stdrun.hoc")
            h.load_file(1, model_name) # http://www.neuron.yale.edu/neuron/static/new_doc/programming/dynamiccode.html#
            h.init()
            h.tstop = tstop
            self.out_data = {}
            self.neurons_names = []
            self.neurons = {}
            #self.sections = {}
            self.__find_all_neurons()
            if len(self.neurons_names) == 0:
                raise RuntimeError(u"In File: {0:s} with model no any neurons has been found. Please check the "
                                   u"the file".format(model_name))
            print self.neurons_names
            for name in self.neurons_names: #TODO put check that we haven't added this neuron yet in dictionary neurons
                self.neurons[name] = MyNeuron(name, index=self.neurons_names.index(name))

            # Initialization of segments and data arrays
            for k, val in self.neurons.iteritems():
                val.init_sections(h, paramVec)
                #for sec in val.sections:
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
        else:
            print 'Simulation is finished'
            sys.exit(0)

    def __find_all_neurons(self):
        """
        Search neurons names from hoc segment name
        """
        from neuron import h
        for h_sec in h.allsec():
            section_name = h_sec.name()
            index = section_name.find('_')
            if index != -1:
                self.neurons_names.append(section_name[0:index])


    def get_time(self):
        from neuron import h
        return h.t

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

    def add_stim(self, n_name, sec_name=10):
        #TODO describe all in detail
        from neuron import h
        stim = h.IClamp(0.5, self.neurons[n_name].sections[sec_name].h_sec)
        stim.amp = 10.0
        stim.delay = 5.0
        stim.dur = 1.0
    def add_synaps(self, sec_id1, sec_id2):
        pass
    def finish(self):
        """
        Do nothing yet
        """
        pass
