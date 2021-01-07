function [eps_BP, sigma_BP] = func_BP_CST_dynamic(freq, target_filename, Cal1_filename, Cal2_filename, resolution)

%%% This is a new version of BP solver which does not need additional parameters (but we do need
%%% three input arguments, the "frequency" (not freq_num now, so you can choose any frequency),
%%% the "target_filename" which heads to the S-parameter file, and "resolution" that will be passed to DOI class

%%% --------------------------------------- Define some basic EM constants ---------------------------------------------

w = 2 * pi * freq;
eps_o = 8.854187817e-12;
uo = 4e-7 * pi;  
source_num = 16;

%%% ----------------- We do the curve fitting here to get EPs of coupling medium and calibration phantoms under arbitary frequency -----------------

freq_array = [700 : 0.5 : 1300] .* 1e6;   %%% These are the inquiry points
fit_freq = [700 : 50 : 1300] .* 1e6;      %%% These are the fitting points
fit_w = 2 * pi * fit_freq;

%%% These are the measured EPs of coupling medium 
eps_r_b = [52.4; 52.06; 52.665; 52.203; 52.1464; 52.1354; 52.8048; 52.6253; 52.0931; 52.4981; 52.102; 51.9399; 51.7879];
sigma_b = [2.9928; 2.5472; 2.7799; 3.3259; 2.9245; 3.2465; 3.0125; 3.9947; 3.859; 3.9225; 4.3962; 4.1691; 4.6537] .* w .* eps_o;

%%% These are the measured EPs of Cal 1 
eps_r_Cal1 = [41.928; 41.2353; 40.8917; 39.45; 39.3708; 38.5727; 37.6621; 37.7461; 36.764; 37.6998; 36.5836; 36.7354; 36.0939];
sigma_Cal1 = [24.8474; 22.0824; 22.9255; 22.0985; 19.9184; 20.525; 19.1796; 19.6627; 17.9339; 18.3372; 17.9049; 16.4417; 16.9719] .* fit_w.' .* eps_o;

%%% These are the measured EPs of Cal 2
eps_r_Cal2 = [54.98; 53.391; 52.2804; 50.829; 50.2496; 49.5723; 48.0418; 47.7047; 46.6832; 46.427; 45.7603; 45.4439; 44.449];
sigma_Cal2 = [25.7777; 24.2249; 24.0116; 23.8182; 22.2428; 23.1322; 22.4505; 22.5571; 21.6656; 20.9165; 21.1112; 20.2122; 20.387] .* fit_w.' .* eps_o;

%%% Fit the EPs of coupling medium
f_eps = fit(fit_freq.' .* 1e-9, eps_r_b, 'poly2');
f_sigma = fit(fit_freq.' .* 1e-9, sigma_b, 'poly2');

fit_coeff_eps = coeffvalues(f_eps);
fit_coeff_sigma = coeffvalues(f_sigma);

fit_eps_r_b = fit_coeff_eps(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_eps(2) .* (freq_array * 1e-9) + fit_coeff_eps(3);
fit_sigma_b = fit_coeff_sigma(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_sigma(2) .* (freq_array * 1e-9) + fit_coeff_sigma(3);

%%% Fit the EPs of Cal 1
f_Cal1_eps = fit(fit_freq.' .* 1e-9, eps_r_Cal1, 'poly2');
f_Cal1_sigma = fit(fit_freq.' .* 1e-9, sigma_Cal1, 'poly2');

fit_coeff_Cal1_eps = coeffvalues(f_Cal1_eps);
fit_coeff_Cal1_sigma = coeffvalues(f_Cal1_sigma);

fit_Cal1_eps = fit_coeff_Cal1_eps(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_Cal1_eps(2) .* (freq_array * 1e-9) + fit_coeff_Cal1_eps(3);
fit_Cal1_sigma = fit_coeff_Cal1_sigma(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_Cal1_sigma(2) .* (freq_array * 1e-9) + fit_coeff_Cal1_sigma(3);

%%% Fit the EPs of Cal 2
f_Cal2_eps = fit(fit_freq.' .* 1e-9, eps_r_Cal2, 'poly2');
f_Cal2_sigma = fit(fit_freq.' .* 1e-9, sigma_Cal2, 'poly2');

fit_coeff_Cal2_eps = coeffvalues(f_Cal2_eps);
fit_coeff_Cal2_sigma = coeffvalues(f_Cal2_sigma);

fit_Cal2_eps = fit_coeff_Cal2_eps(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_Cal2_eps(2) .* (freq_array * 1e-9) + fit_coeff_Cal2_eps(3);
fit_Cal2_sigma = fit_coeff_Cal2_sigma(1) .* (freq_array * 1e-9) .^ 2 + fit_coeff_Cal2_sigma(2) .* (freq_array * 1e-9) + fit_coeff_Cal2_sigma(3);

%%% Now we need to find the corresponding index of freq and define parameters related to EPs of CM, Cal1 and Cal2

[~, freq_idx] = min(abs((freq_array - freq)));

eps_r_b = fit_eps_r_b(freq_idx);
sigma_b = fit_sigma_b(freq_idx);

eps_b = eps_r_b * eps_o - 1i * sigma_b / w;
kb = w * sqrt(uo * eps_b);

Cal1_eps = fit_Cal1_eps(freq_idx);
Cal1_sigma = fit_Cal1_sigma(freq_idx);

Cal2_eps = fit_Cal2_eps(freq_idx);
Cal2_sigma = fit_Cal2_sigma(freq_idx);

%%% -------------------------- Create the DOI class and define compuation configuration ------------------------

DOI_setting = DOI(resolution);

Probes_Tx_x = DOI_setting.ant_xy(:, 2);
Probes_Tx_y = DOI_setting.ant_xy(:, 1);

Probes_Tx_x = Probes_Tx_x * 1e-3;
Probes_Tx_y = Probes_Tx_y * 1e-3;

Probes_Tx = zeros(source_num, 2);
Probes_Tx(:, 1) = Probes_Tx_x;
Probes_Tx(:, 2) = Probes_Tx_y;

Probes_Rx = Probes_Tx;

src_Tx_N = length(Probes_Tx);
src_Rx_N = length(Probes_Rx);

x_dash = DOI_setting.y_axis;
y_dash = DOI_setting.x_axis;

x_dash = x_dash * 1e-3;
y_dash = y_dash * 1e-3;

Nx = length(x_dash);
Ny = length(y_dash);
total_n = Nx * Ny;

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

Axis_x = Axis_x(:);
Axis_y = Axis_y(:);

%%% ---------------------------- Define the incident field --------------------------- %%%

Ez_inc_MoM = zeros(Nx * Ny, src_Tx_N);
Ez_inc_Rx = zeros(src_Rx_N, src_Rx_N);

for src_dash = 1 : src_Tx_N
 
    Probes_pho = sqrt((Axis_x(:) - Probes_Tx(src_dash, 1)) .^ 2 + (Axis_y(:) - Probes_Tx(src_dash, 2)) .^ 2);
    Ez_inc_MoM(:, src_dash) = -1 * w * uo * 0.25 * besselh(0, 2, kb .* Probes_pho');
    
end

[Ess_MoM_Cal1, Ess_MoM_Cal2, Gezz, Gezz_source] = MoM_GreenFunc(freq, resolution, DOI_setting, Ez_inc_MoM, eps_r_b, sigma_b, Cal1_eps, Cal1_sigma, Cal2_eps, Cal2_sigma, src_Tx_N);


%%% ---------------------------- Do the calibration --------------------------- %%%

CST_Cal1_struct = sparameters(Cal1_filename);
CST_Cal2_struct = sparameters(Cal2_filename);

CST_Cal1 = CST_Cal1_struct.Parameters;
CST_Cal2 = CST_Cal2_struct.Parameters;

freq_array_CST = CST_Cal1_struct.Frequencies;

[~, freq_idx] = min(abs((freq_array_CST - freq)));     %%% Find the corresponding frequency index in the CST or measurement data

CST_Cal1 = CST_Cal1(:, :, freq_idx);
CST_Cal2 = CST_Cal2(:, :, freq_idx);

alpha_Ez = (Ess_MoM_Cal1 - Ess_MoM_Cal2) ./ (CST_Cal1 - CST_Cal2);
belta_Ez = Ess_MoM_Cal1 - alpha_Ez .* CST_Cal1;

CST_case_struct = sparameters(target_filename);
CST_case = CST_case_struct.Parameters;

S_scat = alpha_Ez .* CST_case + belta_Ez;


%%% ---------------------------- The BP solver --------------------------- %%%

CST_model_mask = DOI_setting.mask;
zero_index = find(CST_model_mask == 0);

wr_bp_0 = zeros(total_n, src_Tx_N);
u_0 = zeros(total_n, src_Tx_N);
Ge_star_fi = zeros(total_n, src_Tx_N);
Ge_Ge_star_fi = zeros(src_Rx_N, src_Tx_N);

for src = 1 : src_Tx_N
    Ge_star_fi(:, src) = Gezz_source' * S_scat(:, src);
    Ge_Ge_star_fi(:, src) = Gezz_source * Ge_star_fi(:, src);

    wr_bp_0(:, src) = (norm(Ge_star_fi(:, src)) ^ 2 / norm(Ge_Ge_star_fi(:, src)) ^ 2) * Ge_star_fi(:, src);

    u_0(:, src) = Ez_inc_MoM(:, src) + Gezz * wr_bp_0(:, src);
end

X_0 = sum(wr_bp_0 .* conj(u_0), 2) ./ sum(abs(u_0) .^ 2, 2);
X_0(zero_index) = 0;
X_rec_new = reshape(X_0, Nx, Ny);

aaa = X_rec_new;
sigma_BP = -1 .* imag((aaa + 1) * eps_b) * w;
eps_BP = real((aaa + 1) * eps_b) ./ eps_o;

sigma_BP(zero_index) = 0;
eps_BP(zero_index) = 0;































