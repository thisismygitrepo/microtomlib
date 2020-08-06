function [eps_BIM, sigma_BIM, error_obj, error_tmp, XX_array] = BIM_newReg_v1(freq_num, exp_target_filename, exp_Cal1_filename, exp_Cal2_filename, database_path, reg_path)

%%% The DOI has the physical range of x: -100mm : 100mm, y: -115mm : 115mm
load(database_path + "/doi_masks.mat")

%%% Load the regulizer
load(reg_path + "/tik_matrix_1500d_4mm_800MHz.mat")

%%% Load parameters
load(database_path + "/Ez_tot_MoM_array_Cal1_84freq_Alex.mat")
load(database_path + "/Ez_tot_MoM_array_Cal2_84freq_Alex.mat")
load(database_path + "/freq_array_CST_Alex.mat")

w = 2 * pi * freq_array(freq_num);
eps_o = 8.854187817e-12;
uo = 4e-7 * pi;


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


%%% ------------------------------ Adjust the regularizer -------------------------- %%%

tik_diag = diag(Tik);
tik_diag = reshape(tik_diag, Ny, Nx);
tik_diag = tik_diag.';
tik_diag = tik_diag(:);
tik_diag_mat = diag(tik_diag);

average_vector = reshape(average_vector, Ny, Nx);
average_vector = average_vector.';
average_vector = average_vector(:);

eps_complex_tmp = tik_diag_mat * average_vector;

eps_tmp = real(eps_complex_tmp);
sigma_tmp = -1 * imag(eps_complex_tmp) * w * eps_o;

zero_index = find(doi_mask_4mm == 0);
eps_tmp(zero_index) = eps_r_b;
sigma_tmp(zero_index) = sigma_b;
eps_tmp = eps_tmp * eps_o;

X_tmp = (eps_tmp ./ eps_b) - 1i * (sigma_tmp ./ (w * eps_b)) - 1;

%%% ---------------------------- Build the Green function ----------------------------- %%%

al = 4e-3;
a = al / (sqrt(pi));

total_n = length(Axis_x);

Gezz = zeros(total_n, total_n);

for m = 1 : total_n
    x = Axis_x(m);
    y = Axis_y(m);

    p = sqrt((x - Axis_x(:)).^2 + (y - Axis_y(:)).^2);
    
     %%% Green function for Ez field  
    Gezz(m, :) = (kb ^ 2) .* ((1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);
    Gezz(m, m) = (kb ^ 2) * ((1 / (kb ^ 2)) + (1i * pi * a * besselh(1, 2, kb * a)) / (2 * kb));
end

Gezz_source = zeros(src_Rx_N, total_n);
for m = 1 : src_Rx_N
    x = Probes_Rx(m, 1);
    y = Probes_Rx(m, 2);

    p = sqrt((x - Axis_x).^2 + (y - Axis_y).^2);       
    
    %%% Green function for Ez field  
    Gezz_source(m, :) = (kb ^ 2) .* ((-1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);   
end

Gezz = -1 .* Gezz;

%%% ------------------------------ Calibration the measured S-parameters ------------------------------- %%%

Exp_Cal1_struct = sparameters(exp_Cal1_filename);
Exp_Cal2_struct = sparameters(exp_Cal2_filename);

Exp_Cal1 = Exp_Cal1_struct.Parameters;
Exp_Cal2 = Exp_Cal2_struct.Parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The 101th sample is for 700MHz, 126th sample is for 750, 151th sample is for 800MHz, and 176th sample is for 850MHz

freq_sample = [101 : 3 : 351];     %%% Change this to [101 : 3 : 351] for CST simulations
% freq_sample = [101 : 2 : 251];       %%% This is for 76 frequency samples

Exp_Cal1 = Exp_Cal1(:, :, freq_sample(freq_num));
Exp_Cal2 = Exp_Cal2(:, :, freq_sample(freq_num));

Et_MoM_Cal1 = Ez_tot_MoM_array_Cal1(:, :, freq_num);
Et_MoM_Cal2 = Ez_tot_MoM_array_Cal2(:, :, freq_num);

Ess_MoM_Cal1 = Et_MoM_Cal1 - Ez_inc_Rx;
Ess_MoM_Cal2 = Et_MoM_Cal2 - Ez_inc_Rx;

% alpha_Ez = (Ess_MoM_Cal1 - Ess_MoM_Cal2) ./ (Exp_Cal1 - Exp_Cal2);
% belta_Ez = Ess_MoM_Cal1 - alpha_Ez .* Exp_Cal1;
% 
Exp_case = importdata(exp_target_filename);
Exp_case = Exp_case(:, :, freq_num);
% % Exp_case = Exp_case - Exp_Cal1;
% 
% Ess = alpha_Ez .* Exp_case + belta_Ez;

alpha_Ez = Ess_MoM_Cal2 ./ Exp_Cal2;
Ess = alpha_Ez .* Exp_case;

%%% ----------------------------------------------------- The BIM inverse solver -----------------------------------------------------

XX = 0 * ones(total_n, 1);
XX(zero_index) = 0;
% XX = X_tmp;
Ez_tot_BIM = zeros(total_n, src_Tx_N);
% Lambda = 3e2;
Lambda = 5e2;

error_obj = zeros(20, 1);
error_tmp = zeros(20, 1);

XX_array = zeros(Nx * Ny, 20);

for nn = 1 : 20
%     tic
%     fprintf('The %ith iteration\n', nn);
    
    B = [];
    I = eye(total_n);
    
    GG_Ez_BIM = I + Gezz .* repmat(XX.', total_n, 1);
    
    inv_GG_Ez_BIM = GG_Ez_BIM \ I;
    
    for src_num = 1 : src_Tx_N
        Ez_tot_BIM(:, src_num) = inv_GG_Ez_BIM * Ez_inc_MoM(:, src_num);
%         fprintf('Finished the %ith source...\n', src_num);
    end
    
    for src_num = 1 : src_Rx_N
        B_dash = Gezz_source .* repmat(Ez_tot_BIM(:, src_num).', src_Rx_N, 1);
        B = [B ; B_dash];
    end
    
    Breg = [B ; Lambda * eye(total_n)];
    
    x0 = Lambda .* X_tmp;
%     x0 = Lambda * zeros(total_n, 1);
    Ereg = [Ess(:) ; x0];
    
%     XX = lsqr(Breg, Ereg, 1e-6, 5000);
%     XX = lsmr(B, Ess(:), Lambda);
    XX = lsmr(Breg, Ereg);
    
    XX(zero_index) = 0;
    
%     error_obj(nn) = norm(B * XX - Ess(:)) ^ 2 / norm(Ess(:)) ^ 2;
    error_obj(nn) = sqrt(sum(abs(B * XX - Ess(:)) .^ 2) ./ sum(abs(Ess(:)) .^ 2));
%     error_tmp(nn) = norm(XX - X_tmp) ^ 2 / norm(X_tmp) ^ 2;
    XX_rmse = XX;
    XX_rmse(zero_index) = 1;
    XX_tmp_rmse = X_tmp;
    XX_tmp_rmse(zero_index) = 1;
    error_tmp(nn) = sqrt(sum(abs(XX_rmse - XX_tmp_rmse) .^ 2) ./ sum(abs(XX_tmp_rmse) .^ 2));
    
%     fprintf('...Error: %i\n', error_obj(nn));
    
    XX_rec = reshape(XX, Nx, Ny);
    aaa = imresize(XX_rec, 2);
    
    XX_array(:, nn) = XX_rec(:);
    
    eps_BIM = real((aaa + 1) * eps_b) ./ eps_o;
    sigma_BIM = -1 .* imag((aaa + 1) * eps_b) * w;
%     imagesc(eps_BIM);axis image;axis off;colormap(jet);caxis([25 65]);
%     imagesc(sigma_BIM);axis image;axis off;colormap(jet);caxis([0 0.5]);
%     pause(0.001);
%     toc
end

fprintf('The %ith frequency sample is finished...\n', freq_num)
    
