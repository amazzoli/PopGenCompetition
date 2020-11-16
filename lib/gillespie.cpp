#include "gillespie.h"


GillespieBD::GillespieBD(std::mt19937& generator) :
StocProc{generator} { };


void GillespieBD::run(param& params){

    state = vecd(state_dim());
    weights = vecd(state_dim() * 2);
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
        traj_step = params.d.at("traj_step");

    } catch (std::exception) { throw std::runtime_error("Algorithm parameters not found"); }

    run(init_state, end_condition, traj_step);
}


void GillespieBD::run(vecd& init_state, endc_f& end_condition, int traj_step) {

    int step = 0;
    time = 0;
    state_traj = vec2d(0);
    time_traj = vecd(0);
    if (state.size() != init_state.size())
        throw std::runtime_error("Wrong size of the initial state");
    state = init_state;
    
    while (!end_condition(step, state)) {

        update_weights(weights);
        w_tot = 0;
        for (int i=0; i<weights.size(); i++) w_tot += weights[i];
        if (w_tot == 0)
            throw std::runtime_error("0 total weight, all states are zero");
        std::exponential_distribution<double> exp_dist(w_tot);
        double dt = exp_dist(generator);
        time += dt;

        std::discrete_distribution<int> state_sample_dist (weights.begin(), weights.end());
        int state_sample = state_sample_dist(generator);
        // Birth event
        if (state_sample < state_dim()) 
            state[state_sample]++;
        // Death event
        else 
            //state[state_sample - state_dim()] = std::max(0.0, state[state_sample - state_dim()]-1);
            state[state_sample - state_dim()]--;

        if (step % traj_step == 0){
            state_traj.push_back(state);
            time_traj.push_back(time);
        }
        step++;
    }
}


GillespieLV2::GillespieLV2(const param& params, std::mt19937& generator):
GillespieBD{generator} {
    try{
        rhos[0] = params.vecd.at("rhos")[0];
        rhos[1] = params.vecd.at("rhos")[1];
        fs[0] = params.vecd.at("fs")[0];
        fs[1] = params.vecd.at("fs")[1];
        chis[0] = params.vecd.at("chis")[0];
        chis[1] = params.vecd.at("chis")[1];
        M = params.d.at("M");
    }
    catch (std::exception){
        throw std::runtime_error("Lotka Volterra parameters not found");
    }
};


void GillespieLV2::update_weights(vecd& weights) {
    weights[0] = (fs[0] * rhos[0])*state[0];
    weights[1] = (fs[1] * rhos[1])*state[1];
    double death_coef = 0;
    for (int i=0; i<2; i++) death_coef += state[i] * chis[i] * fs[i];
    weights[2] = (rhos[0] * death_coef / M)*state[0];
    weights[3] = (rhos[1] * death_coef / M)*state[1];
}


GillespiePlot2::GillespiePlot2(const param& params, std::mt19937& generator):
GillespieBD{generator} {

    try{
        betas[0] = params.vecd.at("betas")[0];
        betas[1] = params.vecd.at("betas")[1];
        alpha = params.d.at("alpha");
        M = params.d.at("M");
    }
    catch (std::exception){
        throw std::runtime_error("Plotkin parameters not found");
    }
};


void GillespiePlot2::update_weights(vecd& weights) {
    weights[0] = betas[0]*state[0];
    weights[1] = betas[1]*state[1];
    double N = 0;
    for (int i=0; i<2; i++) N += state[i];
    weights[2] = betas[0] * alpha * ( 1 + N / M )*state[0];
    weights[3] = betas[1] * alpha * ( 1 + N / M )*state[1];
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

endc_a_f endc_a_time(int fin_step){
    endc_a_f endc = endc_a_f{
        [fin_step](const int& t, const int s[], const int& l) {
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

endc_a_f endc_a_time_fixation(int fin_step){
    endc_a_f endc = endc_a_f{
        [fin_step](const int& t, const int states[], const int& l) {
            if (t>fin_step) return true;
            for (int i=0; i<l; i++)
                if (states[i] == 0) return true;
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

endc_a_f endc_a_passage(vecd up_bounds, vecd low_bounds) {
    endc_a_f endc = endc_a_f{
        [up_bounds, low_bounds](const int& t, const int states[], const int& l) {
            for (int s=0; s<l; s++) {
                if (states[s] >= up_bounds[s]) return true;
                if (states[s] <= low_bounds[s]) return true;
            }
            return false;
        }
    };
    return endc;
}


GillespieBD* get_gillespieBD(const param& params, std::mt19937& generator) {

    str alg_name;
    try{ alg_name = params.s.at("process_type"); }
    catch (std::exception){ throw std::runtime_error("Stochastic process type not found"); }

    if (alg_name == "lv2"){
        return new GillespieLV2(params, generator);
    }
    else if (alg_name == "plotkin2"){
		return new GillespiePlot2(params, generator);
    }
    else throw std::invalid_argument( "Invalid stochastic process name" );
}


void GillespieBD::print_traj(str path) const {
    vec2d traj = vec2d(0);
    for (int i=0; i<time_traj.size(); i++) {
        vecd aux_traj = vecd(0);
        aux_traj.push_back(time_traj[i]);
        for (int s=0; s<state_traj[i].size(); s++)
            aux_traj.push_back(state_traj[i][s]);
        traj.push_back(aux_traj);
    }
    vecs labels = vecs{"Time"};
    for (int s=0; s<state_traj[0].size(); s++)
        labels.push_back("s"+std::to_string(s));

    print_2d_traj(traj, labels, path);
}