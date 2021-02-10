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
    state_traj = vec2d {init_state};
    time_traj = vecd {0};
    if (state.size() != init_state.size())
        throw std::runtime_error("Wrong size of the initial state");
    state = init_state;
    
    while (!end_condition(step, state)) {

        update_weights(weights);
        w_tot = 0;
        //for (const double& w : weights) std::cout << w << ",";
        //std::cout << "\n";
        for (const double& w : weights) w_tot += w;
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

        if (traj_step > 0 && step % traj_step == 0){
            state_traj.push_back(state);
            time_traj.push_back(time);
        }
        step++;
    }

    state_traj.push_back(state);
    time_traj.push_back(time);
}


GillespieLV2::GillespieLV2(const param& params, std::mt19937& generator):
GillespieBD{generator} {
    try{
        gr[0] = params.vecd.at("rhos")[0];
        gr[1] = params.vecd.at("rhos")[1];
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
    weights[0] = gr[0]*state[0];
    weights[1] = gr[1]*state[1];
    double death_coef = 0;
    for (int i=0; i<2; i++) death_coef += state[i] * chis[i];
    weights[2] = (gr[0] * death_coef / fs[0] / M)*state[0];
    weights[3] = (gr[1] * death_coef / fs[1] / M)*state[1];
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



GillespieChem2::GillespieChem2(const param& params, std::mt19937& generator):
GillespieBD{generator} {
    try{
        delta[0] = params.vecd.at("deltas")[0];
        delta[1] = params.vecd.at("deltas")[1];
        alpha[0] = params.vecd.at("alphas")[0];
        alpha[1] = params.vecd.at("alphas")[1];
        eta[0] = params.vecd.at("etas")[0];
        eta[1] = params.vecd.at("etas")[1];
        M = params.d.at("M");
    }
    catch (std::exception){
        throw std::runtime_error("Chemostat parameters not found");
    }
};


void GillespieChem2::update_weights(vecd& weights) {
    double R = M/(state[0]*alpha[0] + state[1]*alpha[1]);
    weights[0] = eta[0]*alpha[0]*R*state[0];
    weights[1] = eta[1]*alpha[1]*R*state[1];
    weights[2] = delta[0]*state[0];
    weights[3] = delta[1]*state[1];
}


GillespieChemEvolDelta0::GillespieChemEvolDelta0(const param& params, std::mt19937& generator):
GillespieBD{generator} {
    try{
        delta0 = params.d.at("delta0");

        if (params.vecd.find("alphas") != params.vecd.end())
            alphas = params.vecd.at("alphas");
        else{
            double delta_alpha = params.d.at("delta_alpha");
            double alpha_max = params.d.at("alpha_max");
            alphas = vecd{delta_alpha};
            while (alphas[alphas.size()-1] < alpha_max)
                alphas.push_back(alphas[alphas.size()-1]+delta_alpha);
        }

        mut_rate = params.d.at("mut_rate");
        alpha0 = params.d.at("alpha0");
        M = params.d.at("M");
    }
    catch (std::exception){
        throw std::runtime_error("Chemostat parameters not found");
    }
};


void GillespieChemEvolDelta0::update_weights(vecd& weights) {

    double R = 0;
    for (int i=0; i<state.size(); i++) R += state[i] * alphas[i];
    R /= M;

    // First type (with boundary condition)
    double eta1 = 1.0 / (alphas[0] + alpha0);
    double eta2 = 1.0 / (alphas[1] + alpha0);
    weights[0] = ((1-mut_rate)*eta1*alphas[0]*state[0] + mut_rate*eta2*alphas[1]*state[1]) / R;
    weights[state_dim()] = delta0*state[0];

    // Intermediate types
    for (int i=1; i<state.size()-1; i++) {
        weights[i] = mut_rate*eta1*alphas[i-1]*state[i-1]; 
        weights[i] += (1-2*mut_rate)*eta2*alphas[i]*state[i]; 
        eta1 = eta2;
        eta2 = 1.0 / (alphas[i+1] + alpha0);
        weights[i] += mut_rate*eta2*alphas[i+1]*state[i+1]; 
        weights[i] /= R;

        weights[state_dim()+i] = delta0*state[i];
    }

    // Last type
    weights[state_dim()-1] = mut_rate*eta1*alphas[state_dim()-2]*state[state_dim()-2];
    weights[state_dim()-1] += (1-mut_rate)*eta2*alphas[state_dim()-1]*state[state_dim()-1];
    weights[state_dim()-1] /= R;
    weights[2*state_dim()-1] = delta0*state[state_dim()-1];
}


GillespieChemEvolEta0::GillespieChemEvolEta0(const param& params, std::mt19937& generator):
GillespieBD{generator} {
    try{
        eta0 = params.d.at("eta0");

        if (params.vecd.find("alphas") != params.vecd.end())
            alphas = params.vecd.at("alphas");
        else{
            double delta_alpha = params.d.at("delta_alpha");
            double alpha_max = params.d.at("alpha_max");
            alphas = vecd{delta_alpha};
            while (alphas[alphas.size()-1] < alpha_max)
                alphas.push_back(alphas[alphas.size()-1]+delta_alpha);
        }

        mut_rate = params.d.at("mut_rate");
        alpha0 = params.d.at("alpha0");
        M = params.d.at("M");
    }
    catch (std::exception){
        throw std::runtime_error("Chemostat parameters not found");
    }
};


void GillespieChemEvolEta0::update_weights(vecd& weights) {

    double R = 0;
    for (int i=0; i<state.size(); i++) R += state[i] * alphas[i];
    R /= M;

    // First type (with boundary condition)
    weights[0] = eta0 * ((1-mut_rate)*alphas[0]*state[0] + mut_rate*alphas[1]*state[1]) / R;
    weights[state_dim()] = eta0 * (alpha0+alphas[0]) * state[0];

    // Intermediate types
    for (int i=1; i<state.size()-1; i++) {
        weights[i] = mut_rate*alphas[i-1]*state[i-1]; 
        weights[i] += (1-2*mut_rate)*alphas[i]*state[i]; 
        weights[i] += mut_rate*alphas[i+1]*state[i+1]; 
        weights[i] *= eta0/R;

        weights[state_dim()+i] = eta0*(alpha0+alphas[i])*state[i];
    }

    // Last type
    weights[state_dim()-1] = mut_rate*alphas[state_dim()-2]*state[state_dim()-2];
    weights[state_dim()-1] += (1-mut_rate)*alphas[state_dim()-1]*state[state_dim()-1];
    weights[state_dim()-1] *= eta0/R;
    weights[2*state_dim()-1] = eta0*(alpha0+alphas[state_dim()-1])*state[state_dim()-1];
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
    else if (alg_name == "chem2"){
		return new GillespieChem2(params, generator);
    }
    else if (alg_name == "chem_evold"){
		return new GillespieChemEvolDelta0(params, generator);
    }
    else if (alg_name == "chem_evole"){
		return new GillespieChemEvolEta0(params, generator);
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