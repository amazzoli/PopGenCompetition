#ifndef INV_PROB_H
#define INV_PROB_H


#include "ensemble.h"


struct InvProbInfo {
    veci fix0_pop; // Number of extinctions of the resident population
    vec2d time0_pop; // Times of extinctions of the resident population
    veci fixm_inv; // Number of first passages at k*m of the invaders
    vec2d timem_inv; // Times of first passages at k*m of the invaders
    vec2d time0_inv; // Times of extinctions of the invaders
};


void compute_inv_prob_and_print(SPEnsemble* ensemble, param& params, std::string out_path, std::mt19937& generator);

vec2d generate_init_cond(SPEnsemble* ensemble, param& params, int relax_step, int N_init_cond);

InvProbInfo build_inv_prob_info(SPEnsemble* ensemble, param& params, vec2d& init_cond, const vecd& thresholds, std::mt19937& generator);

double compute_inv_prob(InvProbInfo info, int T);

void print_inv_prob(InvProbInfo info, int T, std::string out_path);


#endif