#include "../lib/inv_prob.h"


int main(int argc, char** argv) {

    if (argc != 3) throw std::runtime_error("Parameter path must be passed during execution");

    // Init random generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    // Parsing the parameters
    std::string path(argv[1]), name(argv[2]), data_dir = "../data/";
    std::string full_dir = data_dir + path + "/" + name;
    param params = parse_param_file(full_dir + "param.txt"); // Def in utils

    // Building the process
    //GillespieBD* alg = gillespie_LV(params, generator);
    StocProc* alg = new GillespieLV2(params, generator);
    SPEnsemble* ensemble = new SPEnsemble(alg, false);

    // Algorithm
    Timer timer;
    std::cout << "\n";
    compute_inv_prob_and_print(ensemble, params, full_dir + "inv_p.txt", generator);

    std::cout << "Switching the populations..\n";
    double aux_f = params.vecd["fs"][0], aux_rho = params.vecd["rhos"][0], aux_chi = params.vecd["chis"][0];
    params.vecd["fs"][0] = params.vecd["fs"][1]; params.vecd["fs"][1] = aux_f;
    params.vecd["rhos"][0] = params.vecd["rhos"][1]; params.vecd["rhos"][1] = aux_rho;
    params.vecd["chis"][0] = params.vecd["chis"][1]; params.vecd["chis"][1] = aux_chi;
    //alg = gillespie_LV(params, generator);
    delete alg;
    alg = new GillespieLV2(params, generator);
    
    ensemble = new SPEnsemble(alg, false);
    
    compute_inv_prob_and_print(ensemble, params, full_dir + "inv_p_switch.txt", generator);
    std::cout << timer.elapsed() << "\n";

    delete alg;
    delete ensemble;

    return 0;
}
