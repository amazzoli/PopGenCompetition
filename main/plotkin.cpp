#include "../lib/ensemble.h"



int main(int argc, char** argv) {

    if (argc != 2)
        throw std::runtime_error("Parameter path must be passed during execution");

    // Init random generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    // Parsing the parameters
    std::string path(argv[1]), data_dir = "../data/";
    param params = parse_param_file(data_dir + path + "/param.txt"); // Def in utils

    StocProc* alg = new GillespiePlot2(params, generator);
    SPEnsemble ensemble = SPEnsemble(alg, params.d.at("N_init_cond"));

     // Observables to compute ober trajectory
    vec_d_vecd_f obs = vec_d_vecd_f {
        // n
        d_vecd_f{ [](const vecd& states) { return states[0]; } },
        // n*n
        d_vecd_f{ [](const vecd& states) { return states[0]*states[0]; } },
        // N
        d_vecd_f{ [](const vecd& states) { return states[0] + states[1]; } },
        // N*N
        d_vecd_f{ [](const vecd& states) { return (states[0] + states[1])*(states[0] + states[1]); } },
        // N*n
        d_vecd_f{ [](const vecd& states) { return states[0]*(states[0] + states[1]); } },
        // x
        d_vecd_f{ [](const vecd& states) { return states[0] / (states[0] + states[1]); } },
        // x*x
        d_vecd_f{ [](const vecd& states) { return states[0]*states[0] / (states[0] + states[1]) / (states[0] + states[1]); } },
        // x / N
        d_vecd_f{ [](const vecd& states) { return states[0] / (states[0] + states[1]) / (states[0] + states[1]); } },
        // x*x / N
        d_vecd_f{ [](const vecd& states) { return states[0]*states[0] / (states[0] + states[1]) / (states[0] + states[1]) / (states[0] + states[1]); } },
    }; 

    ensemble.print_averages(obs, params, data_dir + path + "/obs.txt");
    
    ensemble.print_final_states(params, data_dir + path + "/final_state.txt"); 

    delete alg;

    return 0;
}
