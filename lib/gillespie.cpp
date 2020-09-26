#include "gillespie.h"


GillespieBD::GillespieBD(vecRates birth_rates, vecRates death_rates, std::mt19937& generator):
birth_rates{birth_rates}, death_rates{death_rates}, StocProc{generator} {
    if (birth_rates.size() != death_rates.size())
        throw std::runtime_error("Birth and death rates of different sizes");
    state = vecd(birth_rates.size());
    state_size = state.size();
};


void GillespieBD::run(param& params){
    vecd init_state; 
    endc_f end_condition; 
    double time_scale;
    int traj_step;
    try{

        if (params.s.at("end_cond") == "time")
            end_condition = endc_time(params.d.at("max_steps"));
        else if (params.s.at("end_cond") == "fix&time")
            end_condition = endc_time_fixation(params.d.at("max_steps"));
        else if (params.s.at("end_cond") == "passage")
            end_condition = endc_passage(params.vecd.at("up_bound"), params.vecd.at("low_bound"));
        else
            throw std::runtime_error("End condition not valid");

        init_state = params.vecd.at("init_state");
        time_scale = params.d.at("time_scale");
        traj_step = params.d.at("traj_step");

    } catch (std::exception) { throw std::runtime_error("Algorithm parameters not found"); }

    run(init_state, end_condition, time_scale, traj_step);
}


void GillespieBD::run(vecd init_state, endc_f end_condition, double time_scale, int traj_step) {

    int step = 0;
    time = 0;
    state_traj = vec2d(0);
    time_traj = vecd(0);
    vecd weights = vecd(state_size * 2);
    if (state.size() != init_state.size())
        throw std::runtime_error("Wrong size of the initial state");
    state = init_state;
    
    while (!end_condition(step, state)) {

        double W = 0;
        for(int i=0; i<state_size; i++){
            weights[i] = birth_rates[i](state)*state[i] / time_scale;
            weights[state_size + i] = death_rates[i](state)*state[i] / time_scale;
            W += weights[i] + weights[state_size + i];
        }

        std::exponential_distribution<double> exp_dist(1/W);
        double dt = exp_dist(generator);
        time += dt;

        std::discrete_distribution<int> state_sample_dist (weights.begin(), weights.end());
        int state_sample = state_sample_dist(generator);
        // Birth event
        if (state_sample < state_size) 
            state[state_sample]++;
        // Death event
        else 
            state[state_sample - state_size] = std::max(0.0, state[state_sample - state_size]-1);
        
        if (step % traj_step == 0){
            state_traj.push_back(state);
            time_traj.push_back(time);
        }

        step++;
    }
}


void GillespieBD::print_traj(std::string path) const {

    std::ofstream file_r;
    file_r.open(path);

    file_r << "Time\t";
    for (int s=0; s<state_traj[0].size(); s++)
        file_r << "n" << s << "\t";
    file_r << "\n";

    for (int t=0; t<time_traj.size(); t++){
        file_r << time_traj[t] << "\t";
        for (const int& s: state_traj[t])
            file_r << s << "\t";
        file_r << "\n";
    }
    file_r.close();
}


GillespieBD* gillespie_LV(param& params, std::mt19937& generator){

    vecd rhos, fs, chis;
    double M;
    try{
        rhos = params.vecd.at("rhos");
        fs = params.vecd.at("fs");
        chis = params.vecd.at("chis");
        M = params.d.at("M");
    }
    catch (std::exception){
        throw std::runtime_error("Lotka Volterra parameters not found");
    }

    vecRates birth_rates = vecRates(0);
    vecRates death_rates = vecRates(0);

    for (int s=0; s<rhos.size(); s++){
        d_vecd_f birth_rate = d_vecd_f{  [fs, rhos, s](vecd state) {  return fs[s] * rhos[s]; } };
        d_vecd_f death_rate = d_vecd_f{  [fs, rhos, chis, s, M](vecd state) {  
            double death_coef = 0;
            for (int i=0; i<rhos.size(); i++) death_coef += state[i] * chis[i] * fs[i];
            return rhos[s] * death_coef / M; 
        } };
        birth_rates.push_back(birth_rate);
        death_rates.push_back(death_rate);
    }

    return new GillespieBD(birth_rates, death_rates, generator);
}


endc_f endc_time(int fin_step){
    endc_f endc = endc_f{
        [fin_step](int t, vecd s) {
            if (t>fin_step) return true;
            return false;
        }
    };
    return endc;
}


endc_f endc_time_fixation(int fin_step){
    endc_f endc = endc_f{
        [fin_step](int t, vecd states) {
            if (t>fin_step) return true;
            for (const double& state : states)
                if (state == 0) return true;
            return false;
        }
    };
    return endc;
}


endc_f endc_passage(vecd up_bounds, vecd low_bounds) {
    endc_f endc = endc_f{
        [up_bounds, low_bounds](int t, vecd states) {
            for (int s=0; s<states.size(); s++) {
                if (states[s] >= up_bounds[s]) return true;
                if (states[s] <= low_bounds[s]) return true;
            }
            return false;
        }
    };
    return endc;
}