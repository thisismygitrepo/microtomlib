function [eps, sig] = run_multiFreq_solvers_parallel(solver, file_path, file_path_cal1, file_path_cal2, databse_path)

freq_N = 84; 
eps = zeros(58, 50, freq_N);
sig = zeros(58, 50, freq_N);

if solver == "BIM"
    lambda = (6 : (17 - 6) / 83 : 17) * 1e2;   %%% Adaptive lambda for traditional Tikhonov
    itr_num = 15;                              %%% The iteration number for BIM is usually 15~25
    addpath("BIM_solver")
elseif solver == "DBIM"
    lambda = (20 : (90 - 20) / 83 : 100) * 1e2;   %%% Adaptive lambda for traditional Tikhonov
    itr_num = 15;                              %%% The iteration number for DBIM is usually 15~25
    addpath("DBIM_solver")
elseif solver == "CSI"
    lambda = false;
    itr_num = 300;                             %%% The iteration number for CSI is usually 300~1000
    addpath("CSI_solver")
elseif solver == "MRCSI"
    lambda = false;
    itr_num = 300;                             %%% The iteration number for MR-CSI is usually 300~1000
    addpath("MR_CSI_solver")
end


parpool(21)
parfor kk = 1 : freq_N
    tic
    [eps(:, :, kk), sig(:, :, kk)] = BIM_newReg_v1(kk, file_path, file_path_cal1, file_path_cal2, lambda, databse_path, itr_num);
    toc
end


end