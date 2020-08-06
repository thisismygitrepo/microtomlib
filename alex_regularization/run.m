

clc
clear

machine = "Alex";

if machine == "Alex"
    database_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\Parameters_for_solvers";
    reg_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\regularizers";
    file_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\large_stroke_simulations\example1\wIXI087-Guys-0768-T1_t0_1.s16p";
    cal1_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\cals\Cal1.s16p";
    cal2_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\cals\Cal2.s16p";

elseif machine == "Lei"
    database_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/Parameters_for_solvers";
    reg_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/bp_postprocessing";
    file_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/simulations/wiener/large_stroke_simulations/example2/wIXI073-Guys-0755-T1_t3_1.s16p";
    cal1_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/simulations/wiener/cals/Cal1.s16p";
    cal2_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/simulations/wiener/cals/Cal2.s16p";
end

% before running this, make sure your current directory is @ toml/solvers
[eps, sig] = BIM_newReg_v1(1, file_path, cal1_path, cal2_path, database_path, reg_path);



% addpath('/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/bp_postprocessing')
% addpath('/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/bp_postprocessing/Parameters')
% addpath('/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/regularizers')
