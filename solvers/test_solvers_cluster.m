
clc
clear
% This code is designed to test the 4 solvers by trying them out on a
% single example, save and view the results.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This is the path setting from Alex's side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

database_path = "../Parameters_for_solvers";
file_path = "../simulations/wiener/large_stroke_simulations/example1/wIXI087-Guys-0768-T1_t0_1.s16p";
cal1_path = "../simulations/wiener/cals/Cal1.s16p";
cal2_path = "../simulations/wiener/cals/Cal2.s16p";

solvers = ["MRCSI", "CSI", "BIM", "DBIM"];
for i=1:4
    solver = solvers(i);
    [eps, sig] = run_multiFreq_solvers(solver, file_path, cal1_path, cal2_path, database_path);
    save("../tmp/" + solver + ".mat", "eps", "sig")
end

