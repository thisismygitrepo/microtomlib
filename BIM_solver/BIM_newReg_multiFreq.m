clc
clear

freq_N = 84; 

eps_BIM = zeros(58 * 2, 50 * 2, freq_N);
sigma_BIM = zeros(58 * 2, 50 * 2, freq_N);

cnt = 1;
file_name_Cal1 = 'Cal1_Alexs_Sim.s16p';                 %%% File name for calibration phantom 1
file_name_Cal2 = 'Cal2_Alexs_Sim.s16p';					%%% File name for calibration phantom 2
file_name_target = 'wIXI073-Guys-0755-T1_t3_1.s16p';	%%% File name for the target 

lambda = [6 : (17 - 6) / 83 : 17] * 1e2;   %%% Adaptive lambda for traditional Tikhonov

for kk = 1 : freq_N
    tic
    
    [eps_BIM(:, :, cnt), sigma_BIM(:, :, cnt)] = BIM_newReg_v1(kk, file_name_Cal1, file_name_Cal2, file_name_target, lambda);
    cnt = cnt + 1;
    
    toc
end