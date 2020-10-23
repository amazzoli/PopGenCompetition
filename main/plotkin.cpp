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
    SPEnsemble Ensemble = SPEnsemble(alg);

    // Observables to compute ober trajectory
    vec_d_vecd_f obs = vec_d_vecd_f {
        // Fraction of the first speciec
        d_vecd_f{ [](const vecd& states) { return states[0] / (states[0] + states[1]); } },
        // Number of individuals
        d_vecd_f{ [](const vecd& states) { return states[0] + states[1]; } }
    };

    Ensemble.print_moments(obs, params, vecd{1.0, 2.0}, (int)params.d.at("N_real"), data_dir + path + "\\");

    delete[] alg;

    return 0;
}
