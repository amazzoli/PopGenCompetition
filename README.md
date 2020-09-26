# PopGenCompetition
Inverting selection in population genetics models with competition


To run the simulations:

- Set the parameters in the notebook file inside the folder *main* by filling the `params` dictionary.

- Save the paramters in the data folder using `ut.write_params`.

- Compile the c++ script in the *main* folder: for **gillespie_test**: *g++ -o gillespie_test.exe gillespie_test.cpp ../lib/gillespie.cpp ../lib/ensemble.cpp ../lib/utils.cpp -std=c++17*. For **invasion_prob**: *g++ -o invasion_prob.exe invasion_prob.cpp ../lib/gillespie.cpp ../lib/ensemble.cpp ../lib/utils.cpp ../lib/inv_prob.cpp -std=c++17*.

- Execute the script by passing also all the required strings: for **gillespie_test** you need the path of the data folder, e.g. *.\gillespie_test.exe ..\data\test\LV\*, for **invasion_prob** you need the the data folder and the specific name of the parameter file, e.g. *.\invasion_prob.exe inv_prob\LV 1_*.

- The output is generated in the same folder of the parameter file, and can be read in the notebook.