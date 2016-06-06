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

    def __init__(self, problem, update_time= 0.1, sample_time=0.01):
        self.problem = problem
        self.sample_time = sample_time
        self.update_time = update_time
        PlotLayer.deployer = self

    def run(self, vehState, obsState, goalState):
        # overrule data of vehicle and obstacles
        for vehicle in self.problem.vehicles:
            vehicle._overrule_state(vehState['state'])
            vehicle._overrule_input(vehState['input'])
            vehicle.set_terminal_conditions(goalState)

        dummy_obstacle = {'position': [-100., -100.], 'velocity': [0., 0.], 'angle': 0.}
        for i, obstacle in enumerate(self.problem.environment.obstacles):
            # environment contains dummy obstacles to account for new obstacles in the environment
            # filter these out
            obstacle_observation = obsState[i] if obsState is not None and i < len(obsState) else dummy_obstacle
            print 'overruling obstacle state for: ', obstacle, 'with position: ', obstacle.signals['position'][:, -1], ' with new position: ', obstacle_observation['position']
            obstacle._overrule_position(obstacle_observation['position'])
            obstacle._overrule_orientation(obstacle_observation['angle'])
            obstacle._overrule_velocity(obstacle_observation['velocity'])
        # solve problem
        self.problem.solve(0., self.update_time)  # current_time and update_time are 0.
        # check if problem was feasible
        return_status = self.problem.problem.stats()['return_status']
        if  return_status != 'Solve_Succeeded':
            if return_status == 'Infeasible_Problem_Detected':
                status = 'infeasible'
                return status, {}, {}, 0.
            elif return_status == 'Maximum_CpuTime_Exceeded':
                status = 'too slow'
                return status, {}, {}, 0.
            else: 
                print 'Unhandled error message from solver'
                status = 'error'
                return status, {}, {}, 0.
        else:
            status = 'feasible'

        # update everything
        self.problem.update(0., self.update_time, self.sample_time)        
        # check termination criteria
        stop = self.problem.stop_criterium(0., self.update_time)
        # return trajectories
        traj_state, traj_input = {}, {}
        for vehicle in self.problem.vehicles:
            traj_state[str(vehicle)] = vehicle.trajectories['state']
            traj_input[str(vehicle)] = vehicle.trajectories['input']
        motionTime = self.problem.get_variable('T',solution=True)[0][0]
        return status, traj_state, traj_input, motionTime

    def time2index(self, time):
        Ts = self.sample_time
        for k, t in enumerate(self.time):
            t = np.round(t, 6)
            if (t <= time) and (time < (t+Ts)) and ((time-t) <= (t+Ts-time)):
                return k
