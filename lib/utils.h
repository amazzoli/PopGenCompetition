#ifndef UTILS_H
#define UTILS_H


#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>


using veci = std::vector<int>;
using vec2i = std::vector<veci>;
using vecd = std::vector<double>;
using vec2d = std::vector<vecd>;
using vec3d = std::vector<vec2d>;
using vecs = std::vector<std::string>;
using d_vecd_f = std::function<double(vecd)>;
using endc_f = std::function<bool(double, vecd)>;
using vecRates = std::vector<d_vecd_f>;
using dictd = std::map<std::string, double>;
using dictvecd = std::map<std::string, vecd>;
using dicts = std::map<std::string, std::string>;


/* Parameters. They can be doubles, vector of doubles or strings */
struct param {
	dictd d;
	dictvecd vecd;
	dicts s;
};


// IN-OUT

param parse_param_file(std::string file_path);

void print_2d_traj(vec2d traj, vecs labels, std::string file_path);


/* Class for measuring the time between the reset and en enlapsed call */
class Timer {
	private:
		using clock_t = std::chrono::high_resolution_clock;
		using second_t = std::chrono::duration<double, std::ratio<1> >;
		std::chrono::time_point<clock_t> m_beg;
	public:
		Timer() : m_beg(clock_t::now()) { };
		/* Set the onset for the time measure */
		void reset() { m_beg = clock_t::now(); }
		/* Get the time in seconds enlapsed from reset */
		double elapsed() const { return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count(); }
};


class Perc {
	private:
		int m_perc_step;
		int m_max_steps;
		double m_last_perc;
	public:
		Perc(int perc_step, int max_steps) : m_max_steps{max_steps}, m_last_perc{0.0} {
			if (perc_step<1 || perc_step>100) {
				std::cout << "Invalid percentage step. Set by default to 10%\n";
				m_perc_step = 10;
			}
			else m_perc_step = perc_step;
		};
		void step(int curr_step) {
			double perc = (double)curr_step/(double)m_max_steps*100;
			if (perc >= m_last_perc){
				std::cout << round(perc) << "%\n";
				m_last_perc = round(perc) + m_perc_step;
			}
		}
};

#endif