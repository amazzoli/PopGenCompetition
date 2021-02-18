#include "inv_prob.h"


vec2d generate_init_cond(SPEnsemble* ensemble, param& params, int relax_step, int N_init_cond) {

    params.s["end_cond"] = "time";
    params.d["traj_step"] = relax_step-1;
    params.d["max_steps"] = relax_step;

    // Initial conditions at carrying capacity of the second species
    params.vecd["init_state"] = vecd { 0.0, (double)((int)params.d.at("carrying_cap")) };
    std::tuple<vecd, vec2d> init_cond_t = (*ensemble).get_final_states(params);

    vec2d init_cond = vec2d(N_init_cond, vecd(0));
    for (int r=0; r<N_init_cond; r++) {
        //std::cout << std::get<1>(init_cond_t)[r][1] << "\n";
        init_cond[r] = vecd{ 1.0, std::get<1>(init_cond_t)[r][1] };
    }
    return init_cond;
}


InvProbInfo build_inv_prob_info(SPEnsemble* ensemble, param& params, vec2d& init_cond, const vecd& thresholds, std::mt19937& generator) {

    int T = init_cond.size();
    InvProbInfo info = { veci(0), vec2d(0), veci(0), vec2d(0), vec2d(0) };
    params.s["end_cond"] = "passage"; // Changing end condition to first passage
    params.vecd["low_bound"] = vecd { 0, 0 }; // Lower boundary for the end condition
    params.vecd["up_bound"] = vecd { 1000000, 1000000 }; // Upper boundary for the end condition, changed during loop
    params.d["traj_step"] = 0;
    
    for (int k=0; k<thresholds.size(); k++) {

        // Change the end condition of the process
        params.vecd["up_bound"] = vecd { (double)((int)thresholds[k]) , 1000000 };

        // Run the processes until first passage at one of the boundaries
        std::tuple<vecd, vec2d> out = (*ensemble).get_final_states(params, init_cond);
        
        // Update all the observables
        vec2d aux_init_cond = vec2d(0);

        std::cout << "\n";
        vecd time0p = vecd(0), timemi = vecd(0), time0i = vecd(0);
        for (int i=0; i<T; i++){

            std::cout << "Res: " << std::get<1>(out)[i][1] << " Inv: " << std::get<1>(out)[i][0] << "\n";

            // Extinction of the resident population
            if (std::get<1>(out)[i][1] <= 0) time0p.push_back(std::get<0>(out)[i]);

            // Threshold hit by the invaders
            else if (std::get<1>(out)[i][0] >= thresholds[k]) {
                timemi.push_back(std::get<0>(out)[i]);
                aux_init_cond.push_back( vecd { (double)((int)thresholds[k]), (double)((int)std::get<1>(out)[i][1]) } );
            }

            // Extinction of the invaders
            else time0i.push_back(std::get<0>(out)[i]);
        }

        info.fix0_pop.push_back(time0p.size());
        info.fixm_inv.push_back(timemi.size());
        info.time0_pop.push_back(time0p);
        info.timem_inv.push_back(timemi);
        info.time0_inv.push_back(time0i);

        std::cout << "Cycle " << k << " starts, threshold: " << thresholds[k] << ", res_ext: " << time0p.size() << ", inv_passed: " << timemi.size();

        if (timemi.size() == 0) break;

        // Amplify the next initial conditions by sampling from the final conditions of invaredrs at threshold
        double av_pop2 = 0;
        std::uniform_int_distribution<int> dist(0, aux_init_cond.size()-1);
        for (int k=0; k<T; k++){
            init_cond[k] = aux_init_cond[dist(generator)];
            av_pop2 += init_cond[k][1] / double(T);
        }
        std::cout << ", av_pop2 :"  << av_pop2 << "\n"; 
    }

    return info;
}


double compute_inv_prob(InvProbInfo info, int T) {
    double p_inv = info.fix0_pop[0] / (double)T;
    for (int i=1; i<info.fix0_pop.size(); i++) {
        double p_aux = info.fix0_pop[i] / (double)T;
        for (int l=0; l<i; l++) p_aux *= info.fixm_inv[l] / (double)T;
        p_inv += p_aux;
    }
    return p_inv;
}


void print_inv_prob(InvProbInfo info, int T, std::string out_path) {
    std::ofstream out;
    out.open(out_path);

    double p_inv = compute_inv_prob(info, T);
    out << p_inv << "\n";

    for (int m=0; m<info.time0_inv.size(); m++) {
        out << info.fix0_pop[m] << "\t" << info.fixm_inv[m] << "\n";
        for (int t=0; t<info.time0_pop[m].size(); t++)
            out << info.time0_pop[m][t] << "\t";
        out << "\n";
        for (int t=0; t<info.timem_inv[m].size(); t++)
            out << info.timem_inv[m][t] << "\t";
        out << "\n";
        for (int t=0; t<info.time0_inv[m].size(); t++)
            out << info.time0_inv[m][t] << "\t";
        out << "\n";
    }

    out.close();
}


void compute_inv_prob_and_print(SPEnsemble* ensemble, param& params, std::string out_path, std::mt19937& generator) {
    int T, relax_step;
    vecd thresholds;
    try {
        T = params.d.at("N_init_cond");
        relax_step = params.d.at("relax_step");
        thresholds = params.vecd.at("thresholds");
    } catch(std::exception) { throw std::runtime_error("Invalid parameters"); }

    std::cout << "Generating initial conditions..\n";
    vec2d init_cond = generate_init_cond(ensemble, params, relax_step, T);
    std::cout << "Starting the algorithm..\n";
    InvProbInfo info = build_inv_prob_info(ensemble, params, init_cond, thresholds, generator);
    std::cout << "Printing the output..\n";
    print_inv_prob(info, T, out_path);
}