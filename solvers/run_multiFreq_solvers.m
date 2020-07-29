function [eps, sig] = run_multiFreq_solvers(solver, file_path, file_path_cal1, file_path_cal2)

freq_N = 84; 
eps = zeros(58, 50, freq_N);
sig = zeros(58, 50, freq_N);

if solver == "BIM"
    lambda = (6 : (17 - 6) / 83 : 17) * 1e2;   %%% Adaptive lambda for traditional Tikhonov
elseif solver == "DBIM"
    lambda = (20 : (90 - 20) / 83 : 100) * 1e2;   %%% Adaptive lambda for traditional Tikhonov
elseif solver == "CSI" || solver == "MRCSI"
    lambda = false;
end


cnt = 1;


for kk = 1 : freq_N
    tic
    
    if solver == "BIM"
        [eps(:, :, cnt), sig(:, :, cnt)] = BIM_newReg_v1(kk, file_path, file_path_cal1, file_path_cal2, lambda);
    elseif solver == "CSI"
        [eps(:, :, cnt), sig(:, :, cnt)] = CSI_newReg_v1(kk, file_path, file_path_cal1, file_path_cal2);
    elseif solver == "DBIM"
        [eps(:, :, cnt), sig(:, :, cnt)] = DBIM_newReg_v1(kk, file_path, file_path_cal1, file_path_cal2, lambda);
    elseif solver == "MRSCI"
        [eps(:, :, cnt), sig(:, :, cnt)] = MR_CSI_newReg_v1(kk, file_path, file_path_cal2, file_path_cal2);
    end
    
    cnt = cnt + 1;
    toc
end

