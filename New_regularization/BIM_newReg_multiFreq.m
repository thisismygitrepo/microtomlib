clc
clear

freq_N = 84; 

eps_DBIM = zeros(58 * 2, 50 * 2, freq_N);
sigma_DBIM = zeros(58 * 2, 50 * 2, freq_N);

cnt = 1;
file_name_Cal1 = 'CAL1.s16p';
file_name_Cal2 = 'CAL2.s16p';
file_name_target = 'bleeding_z00.mat';

error_obj = zeros(20, freq_N);
error_tmp = zeros(20, freq_N);

XX_array = zeros(58 * 50, 20, freq_N);

for kk = 40 : 40
    tic
    
    [eps_DBIM(:, :, cnt), sigma_DBIM(:, :, cnt), error_obj(:, cnt), error_tmp(:, cnt), XX_array(:, :, cnt)] = BIM_newReg_v1(kk, file_name_Cal1, file_name_Cal2, file_name_target);
    cnt = cnt + 1;

    
%     imagesc(eps_DBIM(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([35 65]);
    imagesc(sigma_DBIM(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([0 0.6]);
    pause(0.001);
    
    toc
end