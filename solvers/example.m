

% On Alex Machine these paths work If you have different machine, commet
% this out (Don't delete) then add your own paths.
database_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\Parameters_for_solvers";
file_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\large_stroke_simulations\example1\wIXI087-Guys-0768-T1_t0_1.s16p";
cal1_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\cals\Cal1.s16p";
cal2_path = "C:\Users\s4551072\OneDrive - The University of Queensland\Pycharm\data\simulations\wiener\cals\Cal2.s16p";

eps, sig = run_multiFreq_solvers("DBIM", file_path, cal1_path, cal2_path, database_path);
