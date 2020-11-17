#include "../lib/inv_prob.h"




int main(int argc, char** argv) {

    if (argc != 2) throw std::runtime_error("Parameter path must be passed during execution");

    // Init random generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    std::string path(argv[1]), data_dir = "../data/";
    std::string full_dir = data_dir + path + "/";
    param params1 = parse_param_file(full_dir + "1_param.txt");
    
    int n_processes = params1.d.at("N_processes");

    GillespieBD* alg;
    SPEnsemble* ensemble;
    Timer time;
    for (int pr=1; pr<=n_processes; pr++) {
        
        // Parameters
        param params = parse_param_file(full_dir + std::to_string(pr) + "_param.txt");

        // Building the precess
        alg = get_gillespieBD(params, generator);
        ensemble = new SPEnsemble(alg, params.d.at("N_init_cond"), false);

        // Compute invasion prob
        compute_inv_prob_and_print(ensemble, params, full_dir+ std::to_string(pr) + "_inv_p.txt", generator);
        std::cout << std::to_string(pr) << "/" << std::to_string(n_processes) << " normal completed\n";

        // Switched parameters
        params = parse_param_file(full_dir + std::to_string(pr) + "sw_param.txt");

        // Building the precess
        alg = get_gillespieBD(params, generator);
        ensemble = new SPEnsemble(alg, params.d.at("N_init_cond"), false);

        // Compute invasion prob
        compute_inv_prob_and_print(ensemble, params, full_dir+ std::to_string(pr) + "sw_inv_p.txt", generator);
        std::cout << std::to_string(pr) << "/" << std::to_string(n_processes) << " switched completed\n"; 
    }
    std::cout << "Time: " << time.elapsed() << "\n";
    
    delete alg;
    delete ensemble;

    return 0;
}
