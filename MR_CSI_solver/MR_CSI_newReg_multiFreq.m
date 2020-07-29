clc
clear

freq_N = 84; 

eps_MR_CSI = zeros(58 * 2, 50 * 2, freq_N);
sigma_MR_CSI = zeros(58 * 2, 50 * 2, freq_N);

cnt = 1;
file_name_Cal1 = 'Cal1_Alexs_Sim.s16p';
file_name_Cal2 = 'Cal2_Alexs_Sim.s16p';
file_name_target = 'wIXI087-Guys-0768-T1_t0_1.s16p';

for kk = 1 : freq_N
    tic
    
    [eps_MR_CSI(:, :, cnt), sigma_MR_CSI(:, :, cnt)] = MR_CSI_newReg_v1(kk, file_name_Cal1, file_name_Cal2, file_name_target);
    cnt = cnt + 1;

    
    imagesc(eps_MR_CSI(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([45 65])
%     imagesc(sigma_MR_CSI(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([0.1 0.3]);view([0 -90])
    pause(0.001);
    
    toc
end