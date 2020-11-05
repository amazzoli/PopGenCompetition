#include "utils.h"


vecd str2vecd(str line, str separator, bool sep_at_end) {
    std::size_t sep_pos = line.find(separator);
    if (sep_pos == str::npos) {
        if (sep_at_end)
            throw std::runtime_error(separator + " separator not found in " + line);
        else {
            vecd v = {std::stod(line.substr(0, sep_pos))};
            return v;
        }
    }

    str elem = line.substr(0, sep_pos);
    vecd v = vecd(0);
    try {
        v.push_back(std::stod(elem));
    }
    catch (std::exception& e){
        v.push_back(0);
    }

    while (true){
        std::size_t next_sep_pos = line.find(separator, sep_pos+1);
        if (sep_at_end && next_sep_pos == str::npos) break;
        str elem = line.substr(sep_pos+1, next_sep_pos-sep_pos);
        try{
            v.push_back(std::stod(elem));
        }
        catch (std::exception& e){
            v.push_back(0);
        }
        if (!sep_at_end && next_sep_pos == str::npos) break;
        sep_pos = next_sep_pos;
    }

    return v;
}


param parse_param_file(str file_path){

    dictd paramd;
    dictvecd paramvecd;
    dicts params;
    
    std::ifstream param_file (file_path);
    if (!param_file.is_open())
        throw std::runtime_error("Error in opening the parameter file at "+file_path);

    str line;
    while ( getline (param_file, line) ) {
        std::size_t tab_pos = line.find("\t");
        str key = line.substr(0,tab_pos);
        str value = line.substr(tab_pos+1, str::npos);

        std::size_t comma_pos = value.find(",");
        if (value.find(",") != str::npos){
            paramvecd[key] = str2vecd(value, ",", true); // Parse a vector
        }
        else{
            try {
                double vald = std::stod(value); // Parse a double
                paramd[key] = vald;
            } catch (std::invalid_argument){
                params[key] = value; // Parse a string if stod gives exception
            }
        }
    }
    param_file.close();

    if (paramd.size() == 0 && paramvecd.size() == 0 && params.size() == 0)
        throw std::runtime_error("Empty parameter file");

    return param{paramd, paramvecd, params};
}


void print_2d_traj(vec2d traj, vecs labels, str file_path) {

    std::ofstream out;
    out.open(file_path);

    for (int k=0; k<labels.size(); k++)
        out << labels[k] << "\t";
    out << "\n";
    
    for (int t=0; t<traj.size(); t++){
        for (int k=0; k<labels.size(); k++)
            out << traj[t][k] << "\t";
        out << "\n";
    }

    out.close();
}


void print_traj(vecd traj, str file_path) {
    std::ofstream out;
    out.open(file_path);
    for (int t=0; t<traj.size(); t++){
        out << traj[t] << "\n";
    }
    out.close();
}