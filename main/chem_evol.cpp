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

    GillespieBD* alg = get_gillespieBD(params, generator);
    for (int i=0; i<params.d.at("N_realizations"); i++){
        (*alg).run(params);
        (*alg).print_traj(data_dir + path+"/traj"+std::to_string(i)+".txt");
        std::cout << i+1 << "/" << params.d.at("N_realizations") << std::endl;
    }
    //SPEnsemble ensemble = SPEnsemble(alg, params.d.at("N_realizations"));


     // Observables to compute ober trajectory

    //ensemble.print_averages(obs, params, data_dir + path + "/obs.txt");

    return 0;
}
