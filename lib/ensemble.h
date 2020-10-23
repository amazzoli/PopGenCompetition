#ifndef ENSEMBLE_H
#define ENSEMBLE_H


#include "gillespie.h"
#include <tuple> 


class SPEnsemble {
    private:
         StocProc* process;
         bool verbose;

         void print_output(std::tuple<vec2d, vec3d> output, vecd moments, std::string dir);

    public:
        SPEnsemble(StocProc* process, bool verbose=true);
         ~SPEnsemble() { delete[] process; }

        /* Get the moments specified in "moments". The parameters should be able to initialize the run 
           method of the StocProc. It returns the time average and time std trajectory as first tuple
           element and the moment trajectories as the second one. */
        std::tuple<vec2d, vec3d> get_moments(param& params, vecd moments, int N_realizations);
        /* Print the moments on an external file */   
        void print_moments(param& params, vecd moments, int N_realizations, std::string dir);

        /* Get the moments specified in "moments" for some fuctions of the state. The parameters should be able to initialize the run 
           method of the StocProc. It returns the time average and time std trajectory as first tuple
           element and the moment trajectories as the second one. */
        std::tuple<vec2d, vec3d> get_moments(vec_d_vecd_f& functions, param& params, vecd moments, int N_realizations);
        /* Print the moments on an external file */   
        void print_moments(vec_d_vecd_f& functions, param& params, vecd moments, int N_realizations, std::string dir);

         /* Get all the states at the end of each trajectory. It returns the times as first element
           of the tuple and the states as the second one. */
        std::tuple<vecd, vec2d> get_final_states(param& params, int N_realizations);
        /* Get all the states at the end of each trajectory. The initial conditions contained in the parameters
           are overwritten by the specified init_cond. It returns the times as first element
           of the tuple and the states as the second one. */
        std::tuple<vecd, vec2d> get_final_states(param& params, vec2d init_cond);
};


#endif