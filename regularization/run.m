
clc
clear

pm = path_manager();
% database_path = fullfile(home, "emvision/Algorithm/toml_data/solvers/Parameters_for_solvers");
file_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/large_stroke_simulations/example1/wIXI087-Guys-0768-T1_t0_1.s16p");
cal1_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/cals/Cal1.s16p");
cal2_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/cals/Cal2.s16p");

[eps, sig] = BIM_alex_reg(700e6, file_path, cal1_path, cal2_path);


% =============================== multifrequency ==========================
% eps_DBIM = zeros(58 * 2, 50 * 2, freq_N);
% sigma_DBIM = zeros(58 * 2, 50 * 2, freq_N);
% error_obj = zeros(20, freq_N);
% error_tmp = zeros(20, freq_N);
% XX_array = zeros(58 * 50, 20, freq_N);
% cnt = 1;
% for kk = 40 : 40
%     tic
%     
%     [eps_DBIM(:, :, cnt), sigma_DBIM(:, :, cnt), error_obj(:, cnt), error_tmp(:, cnt), XX_array(:, :, cnt)] = BIM_newReg_v1(kk, file_name_Cal1, file_name_Cal2, file_name_target);
%     cnt = cnt + 1;
% 
%     
% %     imagesc(eps_DBIM(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([35 65]);
%     imagesc(sigma_DBIM(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([0 0.6]);
%     pause(0.001);
%     
%     toc
% end
