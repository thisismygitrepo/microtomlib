clc
clear

addpath('/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/bp_postprocessing')

%%% The DOI has the physical range of x: -100mm : 100mm, y: -115mm : 115mm
load doi_masks.mat

%%% Load the regulizer
addpath('/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/regularizers')
load tik_matrix_100d_4mm_800MHz.mat

%%% Load parameters
addpath('/Users/uqlguo3/The University of Queensland/Ahmed Al-Saffar - data/Parameters')
load Ez_tot_MoM_array_Cal1_84freq.mat
load Ez_tot_MoM_array_Cal2_84freq.mat
load freq_array_CST.mat

freq_num = 30;

w = 2 * pi * freq_array(freq_num);
eps_o = 8.854187817e-12;
uo = 4e-7 * pi;  

exp_Cal1_filename = 'Cal1.s16p';
exp_Cal2_filename = 'Cal2.s16p';
exp_target_filename = 'Bleeding_z00.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Curve fitting for the permittivity and conductivity of coupling medium

fit_freq = [700 : 50 : 1300] .* 1e6;
fit_w = 2 * pi * fit_freq;
eps_r_b_meas = [52.4; 52.06; 52.665; 52.203; 52.1464; 52.1354; 52.8048; 52.6253; 52.0931; 52.4981; 52.102; 51.9399; 51.7879];
sigma_b_meas = [2.9928; 2.5472; 2.7799; 3.3259; 2.9245; 3.2465; 3.0125; 3.9947; 3.859; 3.9225; 4.3962; 4.1691; 4.6537] .* fit_w.' .* eps_o;

f_eps = fit(fit_freq.' .* 1e-9, eps_r_b_meas, 'poly2');
f_sigma = fit(fit_freq.' .* 1e-9, sigma_b_meas, 'poly2');

fit_coeff_eps = coeffvalues(f_eps);
fit_coeff_sigma = coeffvalues(f_sigma);

fit_eps_r_b = fit_coeff_eps(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_eps(2) .* (freq_array * 1e-9) + fit_coeff_eps(3);
fit_sigma_b = fit_coeff_sigma(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_sigma(2) .* (freq_array * 1e-9) + fit_coeff_sigma(3);

eps_r_b = fit_eps_r_b(freq_num);
sigma_b = fit_sigma_b(freq_num);

eps_b = eps_r_b * eps_o - 1i * sigma_b / w;
kb = w * sqrt(uo * eps_b);

source_num = 16;

%%% --------------------------- Define the position of the waveguide ----------------------------%%%

Probes_Tx_x = [106.6; 83.5; 51.2; 17.4; -16.8; -51.6; -84.2; -106.1; -106.1; -84.2; -51.6; -16.8; 17.4; 51.2; 83.5; 106.6] .* 1e-3;
Probes_Tx_y = [21.1; 55.6; 76.0; 85.4; 86.4; 77.8; 56.5; 21.2; -21.2; -56.5; -77.8; -86.4; -85.4; -76.0; -55.6; -21.1] .* 1e-3;

Probes_Tx = zeros(source_num, 2);
Probes_Tx(:, 1) = Probes_Tx_x;
Probes_Tx(:, 2) = Probes_Tx_y;

Probes_Rx = Probes_Tx;

src_Tx_N = length(Probes_Tx(:, 1));
src_Rx_N = length(Probes_Rx(:, 1));

%%% --------------------------- Define the computation domain ----------------------------- %%%

x_dash = [-115 : 4 : 115] .* 1e-3;
y_dash = [-100 : 4 : 96] .* 1e-3;

Nx = length(x_dash);
Ny = length(y_dash);

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

Axis_x = Axis_x(:);
Axis_y = Axis_y(:);

%%% -------------------------- Find the Probes_Rx mesh index -------------------------- %%%
tmp_x = ones(Nx * Ny, 1);
tmp_y = ones(Nx * Ny, 1);

Probes_Rx_coord = zeros(source_num, 2);
Probes_Rx_index = zeros(source_num, 2);

for kk = 1 : source_num
    dist = sqrt((Axis_x - Probes_Tx_x(kk) .* tmp_x) .^ 2 + (Axis_y - Probes_Tx_y(kk) .* tmp_y) .^ 2);
    [~, min_pos] = min(dist);
    
    [index_x, index_y] = ind2sub([Nx, Ny], min_pos);
    
    Probes_Rx_coord(kk, 1) = Axis_x(min_pos);
    Probes_Rx_coord(kk, 2) = Axis_y(min_pos);
    
    Probes_Rx_index(kk, 1) = index_x;
    Probes_Rx_index(kk, 2) = index_y;
end
    
%%% ---------------------------- Define the incident field --------------------------- %%%

Ez_inc_MoM = zeros(Nx * Ny, src_Tx_N);
Ez_inc_Rx = zeros(src_Rx_N, src_Rx_N);

for src_dash = 1 : src_Tx_N
 
    Probes_pho = sqrt((Axis_x(:) - Probes_Tx(src_dash, 1)) .^ 2 + (Axis_y(:) - Probes_Tx(src_dash, 2)) .^ 2);
    Ez_inc_MoM(:, src_dash) = -1 * w * uo * 0.25 * besselh(0, 2, kb .* Probes_pho');

    Ez_inc_MoM_rec = reshape(Ez_inc_MoM(:, src_dash), Nx, Ny);

    for mm = 1 : src_Rx_N
        Ez_inc_Rx(mm, src_dash) = Ez_inc_MoM_rec(Probes_Rx_index(mm, 1), Probes_Rx_index(mm, 2));
    end
    
end


% %%% ---------------------------- Build the Green function ----------------------------- %%%
% 
% al = 4e-3;
% a = al / (sqrt(pi));
% 
% total_n = length(Axis_x);
% 
% Gezz = zeros(total_n, total_n);
% 
% for m = 1 : total_n
%     x = Axis_x(m);
%     y = Axis_y(m);
% 
%     p = sqrt((x - Axis_x(:)).^2 + (y - Axis_y(:)).^2);
%     
%      %%% Green function for Ez field  
%     Gezz(m, :) = (kb ^ 2) .* ((1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);
%     Gezz(m, m) = (kb ^ 2) * ((1 / (kb ^ 2)) + (1i * pi * a * besselh(1, 2, kb * a)) / (2 * kb));
% end
% 
% Gezz_source = zeros(src_Rx_N, total_n);
% for m = 1 : src_Rx_N
%     x = Probes_Rx(m, 1);
%     y = Probes_Rx(m, 2);
% 
%     p = sqrt((x - Axis_x).^2 + (y - Axis_y).^2);       
%     
%     %%% Green function for Ez field  
%     Gezz_source(m, :) = (kb ^ 2) .* ((-1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);   
% end
% 
% Gezz = -1 .* Gezz;
% 
% %%% ------------------------------ Calibration the measured S-parameters ------------------------------- %%%
% 
% Exp_Cal1_struct = sparameters(exp_Cal1_filename);
% Exp_Cal2_struct = sparameters(exp_Cal2_filename);
% 
% Exp_Cal1 = Exp_Cal1_struct.Parameters;
% Exp_Cal2 = Exp_Cal2_struct.Parameters;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% The 101th sample is for 700MHz, 126th sample is for 750, 151th sample is for 800MHz, and 176th sample is for 850MHz
% 
% freq_sample = [101 : 3 : 351];     %%% Change this to [101 : 3 : 351] for CST simulations
% 
% Exp_Cal1 = Exp_Cal1(:, :, freq_sample(freq_num));
% Exp_Cal2 = Exp_Cal2(:, :, freq_sample(freq_num));
% 
% load Ez_tot_MoM_array_Cal1_84freq.mat     %%% Change this for CST simulations becasue we used 84 frequency samples
% load Ez_tot_MoM_array_Cal2_84freq.mat     %%% Change this for CST simulations becasue we used 84 frequency samples
% 
% Et_MoM_Cal1 = Ez_tot_MoM_array_Cal1(:, :, freq_num);
% Et_MoM_Cal2 = Ez_tot_MoM_array_Cal2(:, :, freq_num);
% 
% Ess_MoM_Cal1 = Et_MoM_Cal1 - Ez_inc_Rx;
% Ess_MoM_Cal2 = Et_MoM_Cal2 - Ez_inc_Rx;
% 
% % alpha_Ez = (Ess_MoM_Cal1 - Ess_MoM_Cal2) ./ (Exp_Cal1 - Exp_Cal2);
% % belta_Ez = Ess_MoM_Cal1 - alpha_Ez .* Exp_Cal1;
% % 
% Exp_case = importdata(exp_target_filename);
% Exp_case = Exp_case(:, :, freq_num);
% % % Exp_case = Exp_case - Exp_Cal1;
% % 
% % Ess = alpha_Ez .* Exp_case + belta_Ez;
% 
% alpha_Ez = Ess_MoM_Cal2 ./ Exp_Cal2;
% Ess = alpha_Ez .* Exp_case;
% 
% %%% ----------------------------------------------------- The DBIM inverse solver -----------------------------------------------------
