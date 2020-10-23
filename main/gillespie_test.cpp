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

    Timer timer;
    StocProc* alg = new GillespieLV2(params, generator);
    (*alg).run(params);
    std::cout << timer.elapsed() << "\n";
    //(*alg).print_traj(data_dir + path + "\\traj.txt");

    timer.reset();
    SPEnsemble LVEnsemble = SPEnsemble(alg);
    LVEnsemble.print_moments(params, vecd{1.0, 2.0}, (int)params.d.at("N_real"), data_dir + path + "\\");
    std::cout << timer.elapsed() << "\n";

    delete alg;

    return 0;
}
