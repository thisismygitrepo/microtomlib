clc
clear

freq_N = 84; 

eps_CSI = zeros(58 * 2, 50 * 2, freq_N);
sigma_CSI = zeros(58 * 2, 50 * 2, freq_N);

cnt = 1;
file_name_Cal1 = 'Cal1_Alexs_Sim.s16p';
file_name_Cal2 = 'Cal2_Alexs_Sim.s16p';
file_name_target = 'wIXI060-Guys-0709-T1_t3_1.s16p';

for kk = 1 : freq_N
    tic
    
    [eps_CSI(:, :, cnt), sigma_CSI(:, :, cnt)] = CSI_newReg_v1(kk, file_name_Cal1, file_name_Cal2, file_name_target);
    cnt = cnt + 1;

    
    imagesc(eps_CSI(:, :, cnt - 1));axis image;axis off;colormap(jet);
%     imagesc(sigma_CSI(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([0.1 0.3]);view([0 -90])
    pause(0.001);
    
    toc
end