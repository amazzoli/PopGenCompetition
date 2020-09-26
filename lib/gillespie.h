#ifndef GILLESPIE_H
#define GILLESPIE_H


#include "utils.h"


/* Abstract class for a stochastic process */
class StocProc {

    protected:
        /* Random number generator */
        std::mt19937 generator;
        /* Trajectory of states generated through run method. First dim: time, second: state space */
        vec2d state_traj;
        /* Trajectory of times of each state generated through run method */
        vecd time_traj;

    public:
        StocProc(std::mt19937& generator) : generator{generator} 
        { state_traj = vec2d(0), time_traj = vecd(0); };

        /* Abstract. Run of the process */
        virtual void run(param& params) = 0;
        /* Abstract. Get dimension of the state space */
        virtual const int state_dim() const = 0;
        /* Get state trajectories */
        const vec2d& get_state_traj() const { return state_traj; }
        /* Get time trajectories */
        const vecd& get_time_traj() const { return time_traj; }
};


/* Gillespie algorithm for birth death processes */
class GillespieBD : public StocProc {

    private:
        /* Per-capita birth rates as list of functions having the state as argument */
        vecRates birth_rates;
        /* Per-capita birth rates as list of functions having the state as argument */
        vecRates death_rates;
        /* Current state */
        vecd state;
        /* Current time */
        double time;
        /* Size of the state space */
        int state_size;

    public:
        GillespieBD(vecRates birth_rates, vecRates death_rates, std::mt19937& generator);
        /* Run give the initial states, end condition, time scale, and number of step after which the stored trajectory is updated*/
        void run(vecd init_state, endc_f end_condition, double time_scale, int traj_step);
        /* Override. Run having all the paramenters sepficied in the structure */
        void run(param& params);
        /* Print the stored state trajectory */
        void print_traj(std::string out_path) const;
        /* Override. Get the size of the state space */
        const int state_dim() const { return state_size; }
};


/* Gillespie for Lotka Volterra */
GillespieBD* gillespie_LV(param& params, std::mt19937& generator);

/* Terminal condition after fin_time steps */
endc_f endc_time(int fin_step);

/* Terminal condition after fin_time steps or if one of the species reaches fixation */
endc_f endc_time_fixation(int fin_step);

/* Terminal condition that stops when a state hits its a bound */
endc_f endc_passage(vecd up_bounds, vecd low_bounds);

#endif