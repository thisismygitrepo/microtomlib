clc
clear


x_dash = [-110 : 2 : 110 + 2] .* 1e-3;  % Hi Lei, what is x_dash mean
y_dash = [-93 : 2 : 92 + 2] .* 1e-3;

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

Probes_Tx_x = [106.6; 83.5; 51.2; 17.4; -16.8; -51.6; -84.2; -106.1; -106.1; -84.2; -51.6; -16.8; 17.4; 51.2; 83.5; 106.6] .* 1e-3;
Probes_Tx_y = [21.1; 55.6; 76.0; 85.4; 86.4; 77.8; 56.5; 21.2; -21.2; -56.5; -77.8; -86.4; -85.4; -76.0; -55.6; -21.1] .* 1e-3;

load eps_bleeding_z00_lambda3_TRA.mat
load sigma_bleeding_z00_lambda3_TRA.mat

% load eps_HH1_CAL2_76freq.mat
% load sigma_HH1_CAL2_76freq.mat

%%% For Patient 1, best freq samples are 1 ~ 30

freq_low = 40;
freq_high = 84;

eps_bleeding = sum(eps_bleeding_z00_lambda3_TRA(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);
sigma_bleeding = sum(sigma_bleeding_z00_lambda3_TRA(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);

% eps_Cal1 = sum(eps_HH1_CAL2_76freq(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);
% sigma_Cal1 = sum(sigma_HH1_CAL2_76freq(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);


% delta_eps_bleed = eps_bleeding - eps_Cal1;
% delta_sigma_bleed = sigma_bleeding - sigma_Cal1;

% delta_sigma_bleed = imgaussfilt(delta_sigma_bleed, 2);
% delta_eps_bleed = imgaussfilt(delta_eps_bleed, 2);

figure
imagesc(imgaussfilt(eps_bleeding, 2));axis image;caxis([45 62]);view([0 -90]);axis off
figure
imagesc(imgaussfilt(sigma_bleeding, 2));axis image;caxis([0 0.6]);view([0 -90]);axis off

% figure
% imagesc(delta_eps_bleed);axis image;view([0 -90]);caxis([-2.5 4]);title('Delta eps bleeding')
% figure
% imagesc(delta_sigma_bleed);axis image;view([0 -90]);caxis([-0.15 0.18]);title('Delta sigma bleeding')
