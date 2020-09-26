#include "ensemble.h"


SPEnsemble::SPEnsemble(StocProc* process, bool verbose) :
process{process}, verbose{verbose} {};


std::tuple<vec2d, vec3d> SPEnsemble::get_moments(param& params, vecd moments, int N_realizations) {

    if (verbose) std::cout << "Computing the moments..\n";

    Perc perc(10, N_realizations);
    vec3d moments_traj = vec3d(0);
    vec2d time_traj = vec2d(2, vecd(0));
    int traj_length = 0;

    for (int r=0; r<N_realizations; r++) {

        if (verbose)
            perc.step(r);

        (*process).run(params);

        if (r==0) {
            // At the first iteration the trajectory length is set ...
            for(double time : (*process).get_time_traj()){
                time_traj[0].push_back(time / (double)N_realizations);
                time_traj[1].push_back(time*time / (double)N_realizations);
            }
            traj_length = time_traj[0].size();
            // ... and the state traj initialized with the proper dimension
            for (int m=0; m<moments.size(); m++)
                moments_traj.push_back( vec2d( traj_length, vecd((*process).state_dim()) ) );
        }
        else {
            // Check on the trajectory length
            vecd sing_time_traj = (*process).get_time_traj();
            if (sing_time_traj.size() != traj_length)
                throw std::runtime_error("Invalid time traj length");
            // Time traj update
            for(int t=0; t<traj_length; t++) {
                time_traj[0][t] += sing_time_traj[t] / (double)N_realizations;
                time_traj[1][t] += sing_time_traj[t]*sing_time_traj[t] / (double)N_realizations;
            }
        }
        // Check on state traj length
        vec2d state_traj = (*process).get_state_traj();
        if (state_traj.size() != traj_length)
                throw std::runtime_error("Invalid state traj length");
        // State traj update
        for (int m=0; m<moments.size(); m++)
            for(int t=0; t<traj_length; t++)
                for (int s=0; s<(*process).state_dim(); s++) {
                    moments_traj[m][t][s] += pow( state_traj[t][s], moments[m] ) / (double)N_realizations;
                }
    }

    if (verbose)
        perc.step(N_realizations);

    return std::tuple<vec2d, vec3d> {time_traj, moments_traj};
}


void SPEnsemble::print_moments(param& params, vecd moments, int N_realizations, std::string dir) {
    std::tuple<vec2d, vec3d> output = get_moments(params, moments, N_realizations);

    for (int m=0; m<moments.size(); m++) {
        char buffer [30];
        std::sprintf(buffer, "moment_%3.2f.txt", moments[m]);

        vecs labels = vecs {"Time_aver", "Time_std"};
        for (int s=0; s<std::get<1>(output)[m][0].size(); s++)
            labels.push_back("x" + std::to_string(s));


        vec2d to_print = vec2d(0);
        for(int t=0; t<std::get<0>(output)[0].size(); t++) {
            double time_av = std::get<0>(output)[0][t];
            double time_m2 = std::get<0>(output)[1][t];
            vecd line = vecd{time_av, sqrt(time_m2 - time_av*time_av)};
            for (int s=0; s<std::get<1>(output)[m][t].size(); s++)
                line.push_back(std::get<1>(output)[m][t][s]);
            to_print.push_back(line);
        }

        print_2d_traj(to_print, labels, dir + buffer);
    }
} 


std::tuple<vecd, vec2d> SPEnsemble::get_final_states(param& params, vec2d init_cond) {

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


std::tuple<vecd, vec2d> SPEnsemble::get_final_states(param& params, int N_realizations) {
    vec2d init_cond = vec2d(N_realizations, params.vecd["init_state"]);
    return get_final_states(params, init_cond);
}
