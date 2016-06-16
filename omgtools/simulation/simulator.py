# This file is part of OMG-tools.
#
# OMG-tools -- Optimal Motion Generation-tools
# Copyright (C) 2016 Ruben Van Parys & Tim Mercy, KU Leuven.
# All rights reserved.
#
# OMG-tools is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

import numpy as np
from plotlayer import PlotLayer


class Simulator:

    def __init__(self, problem, sample_time=0.01, update_time=0.1):
        self.set_problem(problem)
        self.update_time = update_time
        self.sample_time = sample_time
        PlotLayer.simulator = self
        self.reset_timing()

    def set_problem(self, problem):
        self.problem = problem

    def run(self):
        # self.reset_timing()
        self.problem.initialize()
        stop = False
        while not stop:
            stop = self.update()
            if stop:
                update_time = self.problem.vehicles[0].signals['time'][:, -1] - self.current_time
                self.update_timing(update_time)
            else:
                self.update_timing()
        self.problem.final()

    def reset_timing(self):
        self.current_time = 0.
        self.time = np.r_[0.]

    def update_timing(self, update_time=None):
        update_time = self.update_time if not update_time else update_time
        self.current_time += update_time
        n_samp = int(update_time/self.sample_time)
        self.time = np.r_[self.time, np.linspace(
            self.time[-1]+self.sample_time, self.time[-1]+n_samp*self.sample_time, n_samp)]

    def update(self):
        # solve problem
        self.problem.solve(self.current_time, self.update_time)
        # update everything
        self.problem.update(
            self.current_time, self.update_time, self.sample_time)
        # check termination criteria
        stop = self.problem.stop_criterium(self.current_time, self.update_time)
        return stop

    def run_once(self, update=True):
        # initialize problem
        self.problem.initialize()
        # solve problem
        self.problem.solve(0., 0.)
        if not update:
            return None
        # update everything
        self.problem.update(self.current_time, np.inf, self.sample_time)
        self.problem.final()
        # determine timing
        update_time = self.problem.vehicles[0].signals['time'][:, -1] - self.current_time
        self.update_timing(update_time)
        # return trajectories
        trajectories = {}
        for vehicle in self.problem.vehicles:
            trajectories[str(vehicle)] = vehicle.trajectories
        return trajectories

    def time2index(self, time):
        Ts = self.sample_time
        for k, t in enumerate(self.time):
            t = np.round(t, 6)
            if (t <= time) and (time < (t+Ts)) and ((time-t) <= (t+Ts-time)):
                return k
