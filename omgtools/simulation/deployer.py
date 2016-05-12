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


class Deployer:

    def __init__(self, problem, sample_time=0.01):
        self.problem = problem
        self.sample_time = sample_time
        PlotLayer.deployer = self

    def run(self, vehState, obsState, goalState):
        # overrule data of vehicle and obstacles
        for vehicle in self.problem.vehicles:
            vehicle._overrule_state(vehState['state'])
            vehicle._overrule_input(vehState['input'])
            vehicle.set_terminal_conditions(goalState)
        for obstacle in self.problem.environment.obstacles:
            obstacle._overrule_position(obsState[obstacle]['position'])
            obstacle._overrule_orientation(obsState[obstacle]['orientation'])
            obstacle._overrule_velocity(obsState[obstacle]['velocity'])
        # solve problem
        self.problem.solve(0., 0.)  # current_time and update_time are 0.
        # update everything
        self.problem.update(0., 0., self.sample_time)
        # check termination criteria
        stop = self.problem.stop_criterium(0., 0.) #, self.update_time)
        # return trajectories
        traj_state, traj_input = {}, {}
        for vehicle in self.problem.vehicles:
            traj_state[str(vehicle)] = vehicle.trajectories['state']
            traj_input[str(vehicle)] = vehicle.trajectories['input']
        return stop, traj_state, traj_input

    def time2index(self, time):
        Ts = self.sample_time
        for k, t in enumerate(self.time):
            t = np.round(t, 6)
            if (t <= time) and (time < (t+Ts)) and ((time-t) <= (t+Ts-time)):
                return k
