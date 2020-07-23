clc
clear

freq_N = 84; 

eps_DBIM = zeros(58 * 2, 50 * 2, freq_N);
sigma_DBIM = zeros(58 * 2, 50 * 2, freq_N);

cnt = 1;
file_name_Cal1 = 'Cal1_Alexs_Sim.s16p';                        %%% File name for calibration phantom 1
file_name_Cal2 = 'Cal2_Alexs_Sim.s16p';                        %%% File name for calibration phantom 2
file_name_target = 'wIXI060-Guys-0709-T1_t3_1.s16p';           %%% File name for the target

error_obj = zeros(20, freq_N);
error_tmp = zeros(20, freq_N);

XX_array = zeros(58 * 50, 20, freq_N);

lambda = [20 : (90 - 20) / 83 : 100] * 1e2;   %%% Adaptive lambda for traditional Tikhonov

for kk = 1 : freq_N
    tic
    
    [eps_DBIM(:, :, cnt), sigma_DBIM(:, :, cnt), error_obj(:, cnt), error_tmp(:, cnt), XX_array(:, :, cnt)] = DBIM_newReg_v1(kk, file_name_Cal1, file_name_Cal2, file_name_target, lambda);
    cnt = cnt + 1;
    
    toc
end