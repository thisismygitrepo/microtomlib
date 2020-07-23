clc
clear

freq_N = 84; 

eps_DBIM = zeros(58 * 2, 50 * 2, freq_N);
sigma_DBIM = zeros(58 * 2, 50 * 2, freq_N);

cnt = 1;
file_name_Cal1 = 'Cal1_Alexs_Sim.s16p';
file_name_Cal2 = 'Cal2_Alexs_Sim.s16p';
file_name_target = 'wIXI073-Guys-0755-T1_t3_1.s16p';

error_obj = zeros(20, freq_N);
error_tmp = zeros(20, freq_N);

XX_array = zeros(58 * 50, 20, freq_N);

lambda = [6 : (17 - 6) / 83 : 17] * 1e2;   %%% Adaptive lambda for traditional Tikhonov
% lambda = 20e2 * ones(84, 1);     %%% For template regularizer

for kk = 80 : 80
    tic
    
    [eps_DBIM(:, :, cnt), sigma_DBIM(:, :, cnt), error_obj(:, cnt), error_tmp(:, cnt), XX_array(:, :, cnt)] = BIM_newReg_v1(kk, file_name_Cal1, file_name_Cal2, file_name_target, lambda);
    cnt = cnt + 1;

    
%     imagesc(eps_DBIM(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([35 65]);
    imagesc(sigma_DBIM(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([0 1]);view([0 -90])
    pause(0.001);
    
    toc
end