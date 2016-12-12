// This file is part of OMG-tools.

// OMG-tools -- Optimal Motion Generation-tools
// Copyright (C) 2016 Ruben Van Parys & Tim Mercy, KU Leuven.
// All rights reserved.

// OMG-tools is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

#include "Dubins.hpp"
#include <cmath>

using namespace std;

namespace omg{

Dubins::Dubins() : Vehicle(3, 2, 2, 3, 5), poseT(3){

}

void Dubins::setInitialConditions(vector<double>& conditions){
    vector<double> zeros(3);
    setPrediction(conditions, zeros);
    this->pose0 = conditions;
}

void Dubins::setTerminalConditions(vector<double>& conditions){
    this->poseT = conditions;
}

void Dubins::getInitSplineValue(vector<vector<double>>& init_value){
    vector<double> state0(3);
    vector<double> input0(2);
    int n_spl = getNSplines();
    int degree = getDegree();
    int len_basis = getLenBasis();
    getPrediction(state0, input0);

    // Initialize v_til as zeros
    for(int j=0; j<len_basis; j++){
        init_value[0][j] = 0;
    }
    // Initialize tg_ha as linspace between tg_ha0 and tg_haT
    double tg_ha0 = tan(state0[2]/2.);
    double tg_haT = tan(this->poseT[2]/2.);
    for(int j=0; j<len_basis; j++){
        init_value[1][j] = tg_ha0+j*(tg_haT-tg_ha0)/(len_basis-1);
    }
}

void Dubins::setParameters(map<string,vector<double>>& parameters){
    vector<double> state0(3);
    vector<double> input0(2);
    getPrediction(state0, input0);

    vector<double> tg_ha0{tan(state0[2]/2.)};
    parameters["tg_ha0"] = tg_ha0;
    vector<double> v_til0{input0[0]/(1+pow(tg_ha0[0],2))};
    parameters["v_til0"] = v_til0;
    vector<double> dtg_ha0{0.5*input0[1]*(1+pow(tg_ha0[0],2))};
    parameters["dtg_ha0"] = dtg_ha0;
    vector<double> pos0{state0[0],state0[1]};
    parameters["pos0"] = pos0;
    this->pose0[0] = state0[0];
    this->pose0[1] = state0[1];
    this->pose0[2] = state0[2];
    vector<double> posT{this->poseT[0], this->poseT[1]};
    parameters["posT"] = posT;  // x
    vector<double> tg_haT{tan(this->poseT[2]/2.)};
    parameters["tg_haT"] = tg_haT;
}

void Dubins::splines2State(vector<vector<double>>& spline_coeffs, vector<double> time, vector<vector<double>>& state, double sample_time){
    vector<vector<double>> splines(time.size(), vector<double>(spline_coeffs.size())); //size 2 by len_time vector
    this->sampleSplines(spline_coeffs, time, splines);  // uses coeffs to get sampled versions of r and v_til
    vector<double> v_til_s(splines.size());
    vector<double> tg_ha_s(splines.size());
    for (int k=0; k<splines.size(); k++){
        v_til_s[k] = splines[k][0];
        tg_ha_s[k] = splines[k][1];
    }
    vector<double> dx_s(splines.size());
    vector<double> dy_s(splines.size());
    for (int k=0; k<splines.size(); k++){
        dx_s[k] = v_til_s[k]*(1-pow(tg_ha_s[k],2));
        dy_s[k] = v_til_s[k]*(2*tg_ha_s[k]);
    }
    vector<double> int_cte(2); // integration constants
    int_cte[0] = this->pose0[0];//spline_coeffs[0][0];
    int_cte[1] = this->pose0[1];//spline_coeffs[1][0];

    vector<double> x_s(splines.size());
    vector<double> y_s(splines.size());

    this->intNumerically(dx_s, x_s, int_cte[0], sample_time);
    this->intNumerically(dy_s, y_s, int_cte[1], sample_time);

    vector<double> theta_s(splines.size());
    for (int k=0; k<tg_ha_s.size(); k++){
        theta_s[k] = 2*atan2(tg_ha_s[k],1);
    }

    for (int k=0; k<splines.size(); k++){
        state[k][0] = x_s[k];
        state[k][1] = y_s[k];
        state[k][2] = theta_s[k];    
    }    
}

void Dubins::splines2Input(vector<vector<double>>& spline_coeffs, vector<double> time, vector<vector<double>>& input, double sample_time){

    vector<vector<double>> splines(time.size(), vector<double>(spline_coeffs.size()));
    this->sampleSplines(spline_coeffs, time, splines);  // uses coeffs to get sampled versions of r and v_til
    vector<double> v_til_s(splines.size());
    vector<double> tg_ha_s(splines.size());
    for (int k=0; k<splines.size(); k++){
        v_til_s[k] = splines[k][0];
        tg_ha_s[k] = splines[k][1];
    }

    vector<double> dtg_ha_s(splines.size());   
    this->diffNumerically(tg_ha_s, dtg_ha_s, sample_time);
    for (int k=0; k<splines.size(); k++){
        input[k][0] = v_til_s[k]*(1+pow(tg_ha_s[k],2));
        input[k][1] = 2*dtg_ha_s[k]/(1+pow(tg_ha_s[k],2));
    }
}

void Dubins::ode(vector<double>& state, vector<double>& input, vector<double>& dstate){
    // state: x, y, theta
    // inputs: V, dtheta
    // find relation between dstate and state, inputs: dx = Ax+Bu
    // dstate = dx, dy, dtheta
    // dstate[2] = input[1]
    dstate[0] = input[0]*cos(state[2]);
    dstate[1] = input[0]*sin(state[2]);
    dstate[2] = input[1];
}

}