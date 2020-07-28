clc
clear

x_dash = [-115 : 2 : 115 + 2] .* 1e-3;
y_dash = [-100 : 2 : 96 + 2] .* 1e-3;

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

Probes_Tx_x = [106.6; 83.5; 51.2; 17.4; -16.8; -51.6; -84.2; -106.1; -106.1; -84.2; -51.6; -16.8; 17.4; 51.2; 83.5; 106.6] .* 1e-3;
Probes_Tx_y = [21.1; 55.6; 76.0; 85.4; 86.4; 77.8; 56.5; 21.2; -21.2; -56.5; -77.8; -86.4; -85.4; -76.0; -55.6; -21.1] .* 1e-3;

load eps_CSI.mat
load sigma_CSI.mat

% load eps_HH1_CAL2_76freq.mat
% load sigma_HH1_CAL2_76freq.mat

%%% For Patient 1, best freq samples are 1 ~ 30

freq_low = 30;
freq_high = 84;

eps_bleeding = sum(eps_CSI(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);
sigma_bleeding = sum(sigma_CSI(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);

sigma_bleeding = imgaussfilt(sigma_bleeding, 2);

figure
imagesc(imgaussfilt(eps_bleeding, 2));axis image;caxis([45 65]);view([0 -90]);axis off;
figure;imagesc(Axis_y(:), Axis_x(:), sigma_bleeding);axis image;view([0 -90]);axis off;
% hold on; scatter(Probes_Tx_y, Probes_Tx_x, '*b')
% caxis([0.15 0.24])
