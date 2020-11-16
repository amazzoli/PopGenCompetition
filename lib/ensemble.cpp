#include "ensemble.h"


SPEnsemble::SPEnsemble(StocProc* process, int n_realizations, bool verbose) :
process{process}, verbose{verbose}, n_realizations{n_realizations} {};


std::tuple<vec2d, vec2d> SPEnsemble::get_averages(param& params) {
    vec_d_vecd_f funcs = vec_d_vecd_f(0);
    for (int s=0; s<(*process).state_dim(); s++){
        d_vecd_f f = d_vecd_f{ [s](const vecd& states) { return states[s]; } };
        funcs.push_back(f);
    }
    return get_averages(funcs, params);
}


std::tuple<vec2d, vec2d> SPEnsemble::get_averages(vec_d_vecd_f& functions, param& params) {

    if (verbose) std::cout << "Computing the averages..\n";

    Perc perc(10, n_realizations);
    vec2d av_traj = vec2d(0);
    vec2d time_traj = vec2d(2, vecd(0));
    int traj_length = 0;

    for (int r=0; r<n_realizations; r++) {

        if (verbose)
            perc.step(r);

        (*process).run(params);

        if (r==0) {
            // At the first iteration the trajectory length is set ...
            for(double time : (*process).get_time_traj()){
                time_traj[0].push_back(time / (double)n_realizations);
                time_traj[1].push_back(time*time / (double)n_realizations);
            }
            traj_length = time_traj[0].size();
            // ... and the state traj initialized with the proper dimension
            for (int f=0; f<functions.size(); f++)
                av_traj.push_back( vecd(traj_length) );
        }
        else {
            // Check on the trajectory length
            vecd sing_time_traj = (*process).get_time_traj();
            if (sing_time_traj.size() != traj_length)
                throw std::runtime_error("Invalid time traj length");
            // Time traj update
            for(int t=0; t<traj_length; t++) {
                time_traj[0][t] += sing_time_traj[t] / (double)n_realizations;
                time_traj[1][t] += sing_time_traj[t]*sing_time_traj[t] / (double)n_realizations;
            }
        }
        // Check on state traj length
        vec2d state_traj = (*process).get_state_traj();
        if (state_traj.size() != traj_length)
                throw std::runtime_error("Invalid state traj length");
        // State traj update

        for(int t=0; t<traj_length; t++)
            for (int i=0; i<functions.size(); i++) {
                av_traj[i][t] += functions[i](state_traj[t]) / (double)n_realizations;
            }
    }

    if (verbose)
        perc.step(n_realizations);

    return std::tuple<vec2d, vec2d> {time_traj, av_traj};
}


void SPEnsemble::print_averages(vec_d_vecd_f& functions, param& params, str path) {
    std::tuple<vec2d, vec2d> output = get_averages(functions, params);
    print_average_output(output, path);
}


void SPEnsemble::print_averages(param& params, str path) {
    std::tuple<vec2d, vec2d> output = get_averages(params);
    print_average_output(output, path);
} 


std::tuple<vecd, vec2d> SPEnsemble::get_final_states(param& params, vec2d init_cond) {

    if (verbose) std::cout << "Computing final states..\n";
    Perc perc(10, init_cond.size());
    vecd times = vecd(0);
    vec2d states = vec2d(0);

    for (int r=0; r<init_cond.size(); r++) {
        if (verbose) perc.step(r);
        params.vecd["init_state"] = init_cond[r];
        (*process).run(params);
        states.push_back( (*process).get_state_traj()[(*process).get_state_traj().size()-1] );
        times.push_back( (*process).get_time_traj()[(*process).get_time_traj().size()-1] );
    }

    if (verbose) perc.step(init_cond.size());
    return std::tuple<vecd, vec2d> {times, states};
}


std::tuple<vecd, vec2d> SPEnsemble::get_final_states(param& params) {
    vec2d init_cond = vec2d(n_realizations, params.vecd["init_state"]);
    return get_final_states(params, init_cond);
}


void SPEnsemble::print_final_states(param& params, str path){

    std::tuple<vecd, vec2d> output = get_final_states(params);

    vecs labels = vecs {"Time"};
    for (int s=0; s<std::get<1>(output)[0].size(); s++)
        labels.push_back("x" + std::to_string(s));

    vec2d to_print = vec2d(0);
    for(int t=0; t<std::get<0>(output).size(); t++) {
        vecd line = vecd{std::get<0>(output)[t]};
        for (int s=0; s<std::get<1>(output)[t].size(); s++)
            line.push_back(std::get<1>(output)[t][s]);
        to_print.push_back(line);
    }
    
    print_2d_traj(to_print, labels, path);
}


void SPEnsemble::print_average_output(std::tuple<vec2d, vec2d> output, str path){

        vecs labels = vecs {"Time_aver", "Time_std"};
        for (int s=0; s<std::get<1>(output).size(); s++)
            labels.push_back("x" + std::to_string(s));

        vec2d to_print = vec2d(0);
        for(int t=0; t<std::get<0>(output)[0].size(); t++) {
            double time_av = std::get<0>(output)[0][t];
            double time_m2 = std::get<0>(output)[1][t];
            vecd line = vecd{time_av, sqrt(time_m2 - time_av*time_av)};
            for (int s=0; s<std::get<1>(output).size(); s++)
                line.push_back(std::get<1>(output)[s][t]);
            to_print.push_back(line);
        }

        print_2d_traj(to_print, labels, path);
}