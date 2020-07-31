clc
clear

% On Alex Machine these paths work If you have different machine, comment
% this out (Don't delete) then add your own paths.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This is the path setting from Alex's side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% database_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\Parameters_for_solvers";
% file_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\large_stroke_simulations\example1\wIXI087-Guys-0768-T1_t0_1.s16p";
% cal1_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\cals\Cal1.s16p";

% cal2_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\cals\Cal2.s16p";

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This is the path setting from Lei's side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

database_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/Parameters_for_solvers";
file_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/simulations/wiener/large_stroke_simulations/example2/wIXI073-Guys-0755-T1_t3_1.s16p";
cal1_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/simulations/wiener/cals/Cal1.s16p";

cal2_path = "/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/simulations/wiener/cals/Cal2.s16p";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% before running this, make sure your current directory is @ toml/solvers
[eps, sig] = run_multiFreq_solvers("MRCSI", file_path, cal1_path, cal2_path, database_path);
