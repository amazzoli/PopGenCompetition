#ifndef ENSEMBLE_H
#define ENSEMBLE_H


#include "gillespie.h"
#include <tuple> 

using vec_d_vecd_f = std::vector<d_vecd_f>;


/* Stochastic process ensemble. Performs computations over an ensemble of a given stochastic process */
class SPEnsemble {

    private:

         /* Sotchastoc process to analize */
         StocProc* process;
         /* Whether informations are printed on standard output */
         bool verbose;
         /* Number of realizations on the ensemble */
         int n_realizations;

         // AUX METHODS
         void print_average_output(std::tuple<vec2d, vec2d> output, str path);

    public:

        SPEnsemble(StocProc* process, int n_realizations, bool verbose=true);
         ~SPEnsemble() { delete process; }

        /* Get the averages of each dimension. The parameters should be able to initialize the run 
           method of the StocProc. It returns the time average and time std trajectory as first tuple
           element and the averages trajectories as the second one. The stoch proc must generate 
           trajectories of the same lenght, otherwise an error is thown. */
        std::tuple<vec2d, vec2d> get_averages(param& params);
        void print_averages(param& params, str path);

        /* Get the averages for some fuctions of the states. The parameters should be able to initialize the run 
           method of the StocProc. It returns the time average and time std trajectory as first tuple
           element and the average trajectories as the second one.  The stoch proc must generate 
           trajectories of the same lenght, otherwise an error is thown. */
        std::tuple<vec2d, vec2d> get_averages(vec_d_vecd_f& functions, param& params);
        void print_averages(vec_d_vecd_f& functions, param& params, str path);

         /* Get all the states at the end of each trajectory. It returns the times as first element
           of the tuple and the states as the second one. */
        std::tuple<vecd, vec2d> get_final_states(param& params);
        /* Print the final states on an external file */
        void print_final_states(param& params, str path);
        /* Get all the states at the end of each trajectory. The initial conditions contained in the parameters
           are overwritten by the specified init_cond. It returns the times as first element
           of the tuple and the states as the second one. */
        std::tuple<vecd, vec2d> get_final_states(param& params, vec2d init_cond);
};


#endif