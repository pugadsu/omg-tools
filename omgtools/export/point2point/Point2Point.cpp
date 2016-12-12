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

#include "Point2Point.hpp"
#ifdef DEBUG
#include <ctime>
#endif
#include <unistd.h>

using namespace std;
using namespace casadi;

namespace omg{

Point2Point::Point2Point(Vehicle* vehicle,
    double update_time, double sample_time, double horizon_time, int trajectory_length, bool initialize):
parameters(N_PAR), variables(N_VAR), lbg(LBG_DEF), ubg(UBG_DEF) {
    if (trajectory_length > int(horizon_time/sample_time)){
        cerr << "trajectory_length > (horizon_time/sample_time)!" << endl;
    }
    if (trajectory_length > int(update_time/sample_time)){
        time.resize(trajectory_length+1);
        state_trajectory.resize(trajectory_length+1, vector<double>(vehicle->getNState()));
        input_trajectory.resize(trajectory_length+1, vector<double>(vehicle->getNInput()));
    }
    else {
        time.resize(int(update_time/sample_time)+1);
        state_trajectory.resize(int(update_time/sample_time)+1, vector<double>(vehicle->getNState()));
        input_trajectory.resize(int(update_time/sample_time)+1, vector<double>(vehicle->getNInput()));
    }
    this->vehicle = vehicle;
    this->update_time = update_time;
    this->sample_time = sample_time;
    this->horizon_time = horizon_time;
    this->trajectory_length = trajectory_length;
    for (int k=0; k<time.size(); k++){
        time[k] = k*sample_time;
    }
    if (initialize){
        this->initialize();
    }
}

Point2Point::Point2Point(Vehicle* vehicle,
    double update_time, double sample_time, double horizon_time):
Point2Point(vehicle, update_time, sample_time, horizon_time, int(update_time/sample_time), true){
}

Point2Point::Point2Point(Vehicle* vehicle,
    double update_time, double sample_time, double horizon_time, int trajectory_length):
Point2Point(vehicle, update_time, sample_time, horizon_time, trajectory_length, true){
}

void Point2Point::initialize(){
    generateProblem();
    generateSubstituteFunctions();
    args["p"] = parameters;
    cout<<"Initialize"<<endl;
    // cout<<"parameters[0]: "<<parameters[0]<<endl;
    // cout<<"parameters[1]: "<<parameters[1]<<endl;
    // cout<<"parameters[2]: "<<parameters[2]<<endl;
    // cout<<"parameters[3]: "<<parameters[3]<<endl;
    // cout<<"parameters[4]: "<<parameters[4]<<endl;
    // cout<<"parameters[5]: "<<parameters[5]<<endl;
    // cout<<"parameters[6]: "<<parameters[6]<<endl;
    // cout<<"parameters[7]: "<<parameters[7]<<endl;
    // cout<<"parameters[8]: "<<parameters[8]<<endl;
    // cout<<"parameters[9]: "<<parameters[9]<<endl;
    args["x0"] = variables;
    // cout<<"variables[0]: "<<variables[0]<<endl;
    // cout<<"variables[1]: "<<variables[1]<<endl;
    // cout<<"variables[2]: "<<variables[2]<<endl;
    // cout<<"variables[3]: "<<variables[3]<<endl;
    // cout<<"variables[4]: "<<variables[4]<<endl;
    // cout<<"variables[5]: "<<variables[5]<<endl;
    // cout<<"variables[6]: "<<variables[6]<<endl;
    // cout<<"variables[7]: "<<variables[7]<<endl;
    // cout<<"variables[8]: "<<variables[8]<<endl;
    // cout<<"variables[9]: "<<variables[9]<<endl;
    args["lbg"] = lbg;
    args["ubg"] = ubg;
    initSplines();
}

void Point2Point::generateProblem(){
    string obj_path = CASADIOBJ;
    // set options
    Dict options;
    options["ipopt.print_level"] = 0;
    options["print_time"] = 0;
    options["ipopt.tol"] = TOL;
    options["ipopt.linear_solver"] = LINEAR_SOLVER;
    options["ipopt.warm_start_init_point"] = "yes";
    // create nlp solver
    this->problem = nlpsol("problem", "ipopt", obj_path+"/nlp.so", options);
}

void Point2Point::generateSubstituteFunctions(){
@generateSubstituteFunctions@
}

void Point2Point::initSplines(){
@initSplines@
}

void Point2Point::reset(){
    for (int k=0; k<input_trajectory.size(); k++){
        for (int j=0; j<input_trajectory[0].size(); j++){
            input_trajectory[k][j] = 0.0;
        }
    }
}

void Point2Point::resetTime(){
    current_time = 0.0;
    current_time_prev = 0.0;
}

bool Point2Point::update(vector<double>& condition0, vector<double>& conditionT,
    vector<vector<double>>& state_trajectory, vector<vector<double>>& input_trajectory, vector<obstacle_t>& obstacles){
    update(condition0, conditionT, state_trajectory, input_trajectory, obstacles, 0);
}

bool Point2Point::update(vector<double>& condition0, vector<double>& conditionT,
    vector<vector<double>>& state_trajectory, vector<vector<double>>& input_trajectory,
    vector<obstacle_t>& obstacles, int predict_shift){
    #ifdef DEBUG
    double tmeas;
    clock_t begin;
    clock_t end;
    #endif
    // correct current_time with predict_shift:
    current_time += predict_shift*sample_time;
    // transform splines: good init guess for this update
    #ifdef DEBUG
    begin = clock();
    #endif
    transformSplines(current_time, current_time_prev);
    #ifdef DEBUG
    end = clock();
    tmeas = double(end-begin)/CLOCKS_PER_SEC;
    cout << "time in transformSplines: " << tmeas << "s" << endl;
    #endif
    // set target condition
    #ifdef DEBUG
    begin = clock();
    #endif
    vehicle->setTerminalConditions(conditionT);
    #ifdef DEBUG
    end = clock();
    tmeas = double(end-begin)/CLOCKS_PER_SEC;
    cout << "time in setTerminalConditions: " << tmeas << "s" << endl;
    #endif
    // predict initial state and input for problem
    #ifdef DEBUG
    begin = clock();
    #endif
    if (fabs(current_time)<=1.e-6){
        vehicle->setInitialConditions(condition0);
    } else{
        vehicle->predict(condition0, this->state_trajectory, this->input_trajectory, update_time, sample_time, predict_shift);
    }
    #ifdef DEBUG
    end = clock();
    tmeas = double(end-begin)/CLOCKS_PER_SEC;
    cout << "time in predict: " << tmeas << "s" << endl;
    #endif
    // solve problem
    #ifdef DEBUG
    begin = clock();
    #endif
    bool check = solve(current_time, obstacles);
    #ifdef DEBUG
    end = clock();
    tmeas = double(end-begin)/CLOCKS_PER_SEC;
    cout << "time in solve: " << tmeas << "s" << endl;
    #endif
    if (!check){
        current_time_prev = current_time; // prevent to transform again after infeasible!
        return false; // user should retry
    }
    // retrieve splines
    #ifdef DEBUG
    begin = clock();
    #endif
    extractData();
    #ifdef DEBUG
    end = clock();
    tmeas = double(end-begin)/CLOCKS_PER_SEC;
    cout << "time in extractData: " << tmeas << "s" << endl;
    #endif
    // write output state and input trajectory
    for (int k=0; k<trajectory_length; k++){
        for (int j=0; j<state_trajectory[0].size(); j++){
            state_trajectory[k][j] = this->state_trajectory[k][j];
        }
        for (int j=0; j<input_trajectory[0].size(); j++){
            input_trajectory[k][j] = this->input_trajectory[k][j];
        }
    }
    // update current time
    current_time_prev = current_time;
    current_time += update_time;
    return true;
}

bool Point2Point::solve(double current_time, vector<obstacle_t>& obstacles){
    // init variables if first time
    if(fabs(current_time)<=1.e-6){
        initVariables();
    }
    updateBounds(current_time, obstacles);
    setParameters(obstacles);
    args["p"] = parameters;
    cout<<"solve"<<endl;
    // cout<<"parameters[0]: "<<parameters[0]<<endl;
    // cout<<"parameters[1]: "<<parameters[1]<<endl;
    // cout<<"parameters[2]: "<<parameters[2]<<endl;
    // cout<<"parameters[3]: "<<parameters[3]<<endl;
    // cout<<"parameters[4]: "<<parameters[4]<<endl;
    // cout<<"parameters[5]: "<<parameters[5]<<endl;
    // cout<<"parameters[6]: "<<parameters[6]<<endl;
    // cout<<"parameters[7]: "<<parameters[7]<<endl;
    // cout<<"parameters[8]: "<<parameters[8]<<endl;
    // cout<<"parameters[9]: "<<parameters[9]<<endl;

    args["x0"] = variables;
    // cout<<"variables[0]: "<<variables[0]<<endl;
    // cout<<"variables[1]: "<<variables[1]<<endl;
    // cout<<"variables[2]: "<<variables[2]<<endl;
    // cout<<"variables[3]: "<<variables[3]<<endl;
    // cout<<"variables[4]: "<<variables[4]<<endl;
    // cout<<"variables[5]: "<<variables[5]<<endl;
    // cout<<"variables[6]: "<<variables[6]<<endl;
    // cout<<"variables[7]: "<<variables[7]<<endl;
    // cout<<"variables[8]: "<<variables[8]<<endl;
    // cout<<"variables[9]: "<<variables[9]<<endl;
    // cout<<"variables[10]: "<<variables[10]<<endl;
    // cout<<"variables[11]: "<<variables[11]<<endl;
    // cout<<"variables[12]: "<<variables[12]<<endl;
    // cout<<"variables[13]: "<<variables[13]<<endl;
    // cout<<"variables[14]: "<<variables[14]<<endl;
    // cout<<"variables[15]: "<<variables[15]<<endl;
    // cout<<endl;
    // cout<<"variables[16]: "<<variables[16]<<endl;
    // cout<<"variables[17]: "<<variables[17]<<endl;
    // cout<<"variables[18]: "<<variables[18]<<endl;
    // cout<<"variables[19]: "<<variables[19]<<endl;
    // cout<<"variables[20]: "<<variables[20]<<endl;
    // cout<<"variables[21]: "<<variables[21]<<endl;
    // cout<<"variables[22]: "<<variables[22]<<endl;
    // cout<<"variables[23]: "<<variables[23]<<endl;
    // cout<<"variables[24]: "<<variables[24]<<endl;
    // cout<<"variables[25]: "<<variables[25]<<endl;
    // cout<<"variables[26]: "<<variables[26]<<endl;
    // cout<<"variables[27]: "<<variables[27]<<endl;
    // cout<<"variables[28]: "<<variables[28]<<endl;
    // cout<<"variables[29]: "<<variables[29]<<endl;
    // cout<<"variables[30]: "<<variables[30]<<endl;
    // cout<<"variables[31]: "<<variables[31]<<endl;
    // cout<<"variables[32]: "<<variables[32]<<endl;
    // cout<<"variables[33]: "<<variables[33]<<endl;
    // cout<<"variables[34]: "<<variables[34]<<endl;
    // cout<<"variables[35]: "<<variables[35]<<endl;
    // cout<<"variables[36]: "<<variables[36]<<endl;
    // cout<<"variables[37]: "<<variables[37]<<endl;
    // cout<<"variables[38]: "<<variables[38]<<endl;
    // cout<<"variables[39]: "<<variables[39]<<endl;
    // cout<<"variables[40]: "<<variables[40]<<endl;
    // cout<<"variables[41]: "<<variables[41]<<endl;
    // cout<<"variables[41]: "<<variables[42]<<endl;
    // cout<<"variables[41]: "<<variables[43]<<endl;
    // cout<<"variables[41]: "<<variables[44]<<endl;
    // cout<<"variables[41]: "<<variables[45]<<endl;
    // cout<<"variables[41]: "<<variables[46]<<endl;
    // cout<<"variables[41]: "<<variables[47]<<endl;
    // cout<<"variables[41]: "<<variables[48]<<endl;
    // cout<<"variables[41]: "<<variables[49]<<endl;
    // cout<<"variables[41]: "<<variables[50]<<endl;
    // cout<<"variables[41]: "<<variables[51]<<endl;
    // cout<<"variables[41]: "<<variables[52]<<endl;
    // cout<<"variables[41]: "<<variables[53]<<endl;
    // cout<<"variables[41]: "<<variables[54]<<endl;
    // cout<<"variables[41]: "<<variables[55]<<endl;
    // cout<<"variables[41]: "<<variables[56]<<endl;
    // cout<<"variables[41]: "<<variables[57]<<endl;
    // cout<<"variables[41]: "<<variables[58]<<endl;
    // cout<<"variables[41]: "<<variables[59]<<endl;
    // cout<<"variables[41]: "<<variables[60]<<endl;
    // cout<<"variables[41]: "<<variables[61]<<endl;
    // cout<<"variables[41]: "<<variables[62]<<endl;
    // cout<<"variables[41]: "<<variables[63]<<endl;
    // cout<<"variables[41]: "<<variables[64]<<endl;
    // cout<<"variables[41]: "<<variables[65]<<endl;
    // cout<<"variables[41]: "<<variables[66]<<endl;
    // cout<<"variables[41]: "<<variables[67]<<endl;
    // cout<<"variables[41]: "<<variables[68]<<endl;
    // cout<<"variables[41]: "<<variables[69]<<endl;
    // cout<<"variables[41]: "<<variables[70]<<endl;
    // cout<<"variables[41]: "<<variables[71]<<endl;
    // cout<<"variables[41]: "<<variables[72]<<endl;
    // cout<<"variables[41]: "<<variables[73]<<endl;
    // cout<<"variables[41]: "<<variables[74]<<endl;
    // cout<<"variables[41]: "<<variables[75]<<endl;
    // cout<<"variables[41]: "<<variables[76]<<endl;
    // cout<<"variables[41]: "<<variables[77]<<endl;
    // cout<<"variables[41]: "<<variables[78]<<endl;
    // cout<<"variables[41]: "<<variables[79]<<endl;
    // cout<<"variables[41]: "<<variables[80]<<endl;
    // cout<<"variables[41]: "<<variables[81]<<endl;
    // cout<<"variables[41]: "<<variables[82]<<endl;
    // cout<<"variables[41]: "<<variables[83]<<endl;
    // cout<<"variables[41]: "<<variables[84]<<endl;
    // cout<<endl;

    args["lbg"] = lbg;
    // cout<<"lbg[0]: "<<lbg[0]<<endl;
    // cout<<"lbg[1]: "<<lbg[1]<<endl;
    // cout<<"lbg[2]: "<<lbg[2]<<endl;
    // cout<<"lbg[3]: "<<lbg[3]<<endl;
    // cout<<"lbg[4]: "<<lbg[4]<<endl;
    // cout<<"lbg[5]: "<<lbg[5]<<endl;
    // cout<<"lbg[6]: "<<lbg[6]<<endl;
    // cout<<"lbg[7]: "<<lbg[7]<<endl;

    args["ubg"] = ubg;
    // cout<<"ubg[0]: "<<ubg[0]<<endl;
    // cout<<"ubg[1]: "<<ubg[1]<<endl;
    // cout<<"ubg[2]: "<<ubg[2]<<endl;
    // cout<<"ubg[3]: "<<ubg[3]<<endl;
    // cout<<"ubg[4]: "<<ubg[4]<<endl;
    // cout<<"ubg[5]: "<<ubg[5]<<endl;
    // cout<<"ubg[6]: "<<ubg[6]<<endl;
    // cout<<"ubg[7]: "<<ubg[7]<<endl;

    sol = problem(args);
    solver_output = string(problem.stats().at("return_status"));
    if (solver_output.compare("Solve_Succeeded") != 0){
        cout << solver_output << endl;
        return false;
    } else{
        vector<double> var(sol.at("x"));
        for (int k=0; k<n_var; k++){
            variables[k] = var[k];
            cout<<"solved var: "<<variables[k]<<endl;
        }
        return true;
    }
}

void Point2Point::initVariables(){
    map<string, map<string, vector<double>>> var_dict;
    int n_spl = vehicle->getNSplines();
    int len_basis = vehicle->getLenBasis();
    vector<vector<double>> init_var_veh (n_spl, vector<double>(len_basis));
    vehicle->getInitSplineValue(init_var_veh);
    vector<double> init_var_veh_vec(n_spl*len_basis);
    for (int k=0; k<n_spl; k++){
        for (int j=0; j<len_basis; j++){
            init_var_veh_vec[k*len_basis+j] = init_var_veh[k][j];
        }
    }
    var_dict[VEHICLELBL]["splines0"] = init_var_veh_vec;
    getVariableVector(variables, var_dict);
}

void Point2Point::setParameters(vector<obstacle_t>& obstacles){
    map<string, map<string, vector<double>>> par_dict;
    fillParameterDict(obstacles, par_dict);
    getParameterVector(parameters, par_dict);
}

void Point2Point::fillParameterDict(vector<obstacle_t>& obstacles, map<string, map<string, vector<double>>>& par_dict){
    map<string, vector<double>> par_dict_veh;
    vehicle->setParameters(par_dict_veh);
    par_dict[VEHICLELBL] = par_dict_veh;
    if (!freeT){
        par_dict[P2PLBL]["t"] = {fmod(round(current_time*1000.)/1000., horizon_time/(vehicle->getKnotIntervals()))};
        par_dict[P2PLBL]["T"] = {horizon_time};
    } else{
        par_dict[P2PLBL]["t"] = {0.0};
    }
    string obstacle_lbls [N_OBS] = OBSTACLELBLS;
    std::vector<double> pos0(2), vel0(2), acc0(2);
    std::vector<double> posT(2), velT(2), accT(2);
    for (int k=0; k<n_obs; k++){
        pos0 = obstacles[k].position;
        vel0 = obstacles[k].velocity;
        acc0 = obstacles[k].acceleration;
        // prediction over update_time
        for (int j=0; j<2; j++){
            posT[j] = pos0[j] + update_time*vel0[j] + 0.5*pow(update_time,2)*acc0[j];
            velT[j] = vel0[j] + update_time*acc0[j];
            accT[j] = acc0[j];
        }
        par_dict[obstacle_lbls[k]]["x"] = posT;
        par_dict[obstacle_lbls[k]]["v"] = velT;
        par_dict[obstacle_lbls[k]]["a"] = accT;
        par_dict[obstacle_lbls[k]]["checkpoints"] = obstacles[k].checkpoints;
        par_dict[obstacle_lbls[k]]["rad"] = obstacles[k].radii;
    }
}

void Point2Point::extractData(){
    map<string, map<string, vector<double>>> var_dict;
    getVariableDict(variables, var_dict);
    vector<double> spline_coeffs_vec(var_dict[VEHICLELBL]["splines0"]);
    vehicle->setKnotHorizon(horizon_time);
    if (freeT){
        horizon_time = var_dict[P2PLBL]["T"][0];
    }
    vehicle->setKnotHorizon(horizon_time);
    int n_spl = vehicle->getNSplines();
    int len_basis = vehicle->getLenBasis();
    vector<vector<double>> spline_coeffs(n_spl, vector<double>(len_basis));
    for (int k=0; k<n_spl; k++){
        for (int j=0; j<len_basis; j++){
            spline_coeffs[k][j] = spline_coeffs_vec[k*len_basis+j];
        }
    }
    retrieveTrajectories(spline_coeffs);
}

void Point2Point::retrieveTrajectories(vector<vector<double>>& spline_coeffs){
    vector<double> time(this->time);
    if (!freeT){
        for (int k=0; k<time.size(); k++){
            time[k] += fmod(round(current_time*1000.)/1000., horizon_time/vehicle->getKnotIntervals());
        }
    }
    vehicle->splines2State(spline_coeffs, time, state_trajectory, this->sample_time);
    vehicle->splines2Input(spline_coeffs, time, input_trajectory, this->sample_time);
    cout<<"in retrieveTrajectories"<<endl;
}

void Point2Point::getParameterVector(vector<double>& par_vect, map<string, map<string, vector<double>>>& par_dict){
@getParameterVector@
}

void Point2Point::getVariableVector(vector<double>& var_vect, map<string, map<string, vector<double>>>& var_dict){
@getVariableVector@
}

void Point2Point::getVariableDict(vector<double>& var_vect, map<string, map<string, vector<double>>>& var_dict){
@getVariableDict@
}

void Point2Point::updateBounds(double current_time, vector<obstacle_t>& obstacles){
@updateBounds@
}

void Point2Point::transformSplines(double current_time, double current_time_prev){
@transformSplines@
}


}
