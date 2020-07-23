clc
clear

x_dash = [-115 : 2 : 115 + 2] .* 1e-3;
y_dash = [-100 : 2 : 96 + 2] .* 1e-3;

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

Probes_Tx_x = [106.6; 83.5; 51.2; 17.4; -16.8; -51.6; -84.2; -106.1; -106.1; -84.2; -51.6; -16.8; 17.4; 51.2; 83.5; 106.6] .* 1e-3;
Probes_Tx_y = [21.1; 55.6; 76.0; 85.4; 86.4; 77.8; 56.5; 21.2; -21.2; -56.5; -77.8; -86.4; -85.4; -76.0; -55.6; -21.1] .* 1e-3;

load eps_lambda2p5.mat
load sigma_lambda2p5.mat

% load eps_HH1_CAL2_76freq.mat
% load sigma_HH1_CAL2_76freq.mat

%%% For Patient 1, best freq samples are 1 ~ 30

freq_low = 1;
freq_high = 70;

eps_bleeding = sum(eps_DBIM(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);
sigma_bleeding = sum(sigma_DBIM(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);

sigma_bleeding = imgaussfilt(sigma_bleeding, 2);

% min_sigma = min(sigma_bleeding(:));
% max_sigma = max(sigma_bleeding(:));
% 
% sigma_bleeding = (sigma_bleeding - min_sigma) .* (1 / (max_sigma - min_sigma));

% eps_Cal1 = sum(eps_HH1_CAL2_76freq(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);
% sigma_Cal1 = sum(sigma_HH1_CAL2_76freq(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);


% delta_eps_bleed = eps_bleeding - eps_Cal1;
% delta_sigma_bleed = sigma_bleeding - sigma_Cal1;

% delta_sigma_bleed = imgaussfilt(delta_sigma_bleed, 2);
% delta_eps_bleed = imgaussfilt(delta_eps_bleed, 2);

% figure
% imagesc(imgaussfilt(eps_bleeding, 2));axis image;caxis([45 62]);view([0 -90]);axis off;
figure
imagesc(sigma_bleeding);axis image;caxis([0 1]);view([0 -90]);axis off;

% figure
% imagesc(delta_eps_bleed);axis image;view([0 -90]);caxis([-2.5 4]);title('Delta eps bleeding')
% figure
% imagesc(delta_sigma_bleed);axis image;view([0 -90]);caxis([-0.15 0.18]);title('Delta sigma bleeding')


% load error_obj_lambda3_TRA
% 
% freq_num = 50;    %%% select the frequency bin to plot the error curve
% 
% figure;plot(error_obj_lambda3_TRA(:, freq_num));


% [Nx, Ny] = size(eps_DBIM(:, :, 1));
% 
% mask_1 = zeros(Nx, Ny);
% 
% center_x = 0;
% center_y = 0;
% 
% X_radius = 136e-3 / 2;
% Y_radius = 86e-3 / 2;
% 
% for ii = 1 : Nx
%     for jj = 1 : Ny
%         
%         xx = Axis_x(ii, jj);
%         yy = Axis_y(ii, jj);
%         
%         pho = (xx - center_x) ^ 2 / X_radius ^ 2 + (yy - center_y) ^ 2 / Y_radius ^ 2;
%         
%         if pho <= 1
%             mask_1(ii, jj) = 1;
%         else
%             mask_1(ii, jj) = 0;
%         end
%     end
% end
% 
% mask_2 = imresize(mask_1(24 : 93, 30 : 72), 2);
% 
% figure
% imagesc(imresize(imgaussfilt(sigma_bleeding(24 : 93, 30 : 72), 2), 2) .* mask_2);axis image;caxis([0.1 0.25]);axis off;
% view([0 -90]);