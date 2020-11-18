#ifndef GILLESPIE_H
#define GILLESPIE_H


#include "utils.h"


/* Abstract class for a stochastic process */
class StocProc {

    protected:
        /* Random number generator */
        std::mt19937 generator;
        /* Trajectory of states generated through run method. First dim: time step, second: state space */
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


/* Abstract Gillespie algorithm for birth death processes */
class GillespieBD : public StocProc {

    private:

        double w_tot;

    protected:
        /* Current state */
        vecd state;
        /* Current time */
        double time;
        /* Size of the state space */
        //const static int state_d = 2;
        //double weights[2*state_d];
        vecd weights;
        //virtual void update_weights(double (&weights)[2*state_d]) = 0;
        virtual void update_weights(vecd& weights) = 0;
        
    public:
        GillespieBD(std::mt19937& generator);
        /* Run give the initial states, end condition, time scale, and number of step after which the stored trajectory is updated*/
        void run(vecd& init_state, endc_f& end_condition, int traj_step);
        /* Override. Run having all the paramenters sepficied in the structure */
        void run(param& params);
        /* Print the stored state trajectory */
        void print_traj(str out_path) const;
        /* Override. Get the size of the state space */
        virtual const int state_dim() const = 0;
};


/* Gillespie algorithm for 2-species Lotka Volterra  */
class GillespieLV2 : public GillespieBD {

    private:
        /* Pop time scale */
        double rhos[2];
        /* Pop fitness */
        double fs[2];
        /* Pop competitive disadvantage */
        double chis[2];
        /* Pop size */
        int M;

    protected:
    
        void update_weights(vecd& weights);

    public:
        GillespieLV2(const param& params, std::mt19937& generator);
        const int state_dim() const { return 2; }
};


/* Gillespie algorithm for 2-species generalized Moran model (Plotkin 2010) */
class GillespiePlot2 : public GillespieBD {

    private:
        /* Pop time scale */
        double betas[2];
        /* Pop fitness */
        double alpha;
        /* Pop size */
        int M;

    protected:
        void update_weights(vecd& weights);

    public:
        GillespiePlot2(const param& params, std::mt19937& generator);
        const int state_dim() const { return 2; }
};


/* Terminal condition after fin_time steps */
endc_f endc_time(int fin_step);
endc_a_f endc_a_time(int fin_step);

/* Terminal condition after fin_time steps or if one of the species reaches fixation */
endc_f endc_time_fixation(int fin_step);
endc_a_f endc_a_time_fixation(int fin_step);

/* Terminal condition that stops when a state hits its a bound */
endc_f endc_passage(vecd up_bounds, vecd low_bounds);
endc_a_f endc_a_passage(vecd up_bounds, vecd low_bounds);



GillespieBD* get_gillespieBD(const param& params, std::mt19937& generator);



#endif