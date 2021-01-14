function [eps_BP, sigma_BP] = func_BP_CST_dynamic(freq, resolution, target_filename, Cal1_filename, Cal2_filename)


%%% the "target_filename" which heads to the S-parameter file, and "resolution" that will be passed to DOI class

%%% Technical note: Cals are used here for the purpose of matching input
%%% data coming from 3D domain (whether reality or CST simulations) match
%%% the 2D nature of BP algorithm. This is contrast to using them for
%%% matching simulation to reality or the other direction.


%%% --------------------------------------- Define some basic EM constants ---------------------------------------------

pm = path_manager();
const = load(pm.join(pm.gdrive, "emvision\Algorithm\toml_data\system\em_const.mat"));
w = 2 * pi * freq;

%%% ----------------- We do the curve fitting here to get EPs of coupling medium and calibration phantoms under arbitary frequency -----------------

pm = path_manager();
%%% Fit the EPs of coupling medium
[eps_r_b, eps_pp_b] = getProperties(pm.join(pm.gdrive, "/emvision\Algorithm\toml_data\system\Material Properties\cm.txt"), freq);
sigma_b = eps_pp_b * w * const.eps_o;
eps_b = eps_r_b * const.eps_o - 1i * sigma_b / w;
kb = w * sqrt(const.uo * eps_b);

%%% Fit the EPs of Cal 1
[Cal1_eps, Cal1_eps_pp] = getProperties(pm.join(pm.gdrive, "/emvision\Algorithm\toml_data\system\Material Properties\Cal1_Prop.txt"), freq);
Cal1_sigma = Cal1_eps_pp * w * const.eps_o;

%%% Fit the EPs of Cal 2
[Cal2_eps, Cal2_eps_pp] = getProperties(pm.join(pm.gdrive, "/emvision\Algorithm\toml_data\system\Material Properties\Cal2_Prop.txt"), freq);
Cal2_sigma = Cal2_eps_pp * w * const.eps_o;

%%% -------------------------- Create the DOI class and define compuation configuration ------------------------

DOI_setting = DOI(resolution);

Probes_Tx_x = DOI_setting.ant_xy(:, 2);
Probes_Tx_y = DOI_setting.ant_xy(:, 1);

Probes_Tx_x = Probes_Tx_x * 1e-3;
Probes_Tx_y = Probes_Tx_y * 1e-3;

Probes_Tx = zeros(DOI_setting.source_num, 2);
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

for src_dash = 1 : src_Tx_N
 
    Probes_pho = sqrt((Axis_x(:) - Probes_Tx(src_dash, 1)) .^ 2 + (Axis_y(:) - Probes_Tx(src_dash, 2)) .^ 2);
    Ez_inc_MoM(:, src_dash) = -1 * w * const.uo * 0.25 * besselh(0, 2, kb .* Probes_pho');
    
end


green_path = pm.home + "/cache/BPgreen" + string(freq) + ".mat";
green_path = char(green_path);
try
    green = load(green_path);
    green = green.green;
    fprintf('Green function loaded from cache...\n');
catch
    green = MoM_GreenFunc(freq, DOI_setting, Ez_inc_MoM, eps_r_b, sigma_b, Cal1_eps, Cal1_sigma, Cal2_eps, Cal2_sigma, src_Tx_N);
    save(green_path, "green")
end

%%% ---------------------------- Do the calibration --------------------------- %%%

CST_Cal1_struct = sparameters(Cal1_filename);
CST_Cal2_struct = sparameters(Cal2_filename);

CST_Cal1 = CST_Cal1_struct.Parameters;
CST_Cal2 = CST_Cal2_struct.Parameters;

freq_array_CST = CST_Cal1_struct.Frequencies;

[~, freq_idx] = min(abs((freq_array_CST - freq)));     %%% Find the corresponding frequency index in the CST or measurement data

CST_Cal1 = CST_Cal1(:, :, freq_idx);
CST_Cal2 = CST_Cal2(:, :, freq_idx);

alpha_Ez = (green.Ess_Cal1 - green.Ess_Cal2) ./ (CST_Cal1 - CST_Cal2);
belta_Ez = green.Ess_Cal1 - alpha_Ez .* CST_Cal1;

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
    Ge_star_fi(:, src) = green.Gezz_source' * S_scat(:, src);
    Ge_Ge_star_fi(:, src) = green.Gezz_source * Ge_star_fi(:, src);

    wr_bp_0(:, src) = (norm(Ge_star_fi(:, src)) ^ 2 / norm(Ge_Ge_star_fi(:, src)) ^ 2) * Ge_star_fi(:, src);

    u_0(:, src) = Ez_inc_MoM(:, src) + green.Gezz * wr_bp_0(:, src);
end

X_0 = sum(wr_bp_0 .* conj(u_0), 2) ./ sum(abs(u_0) .^ 2, 2);
X_0(zero_index) = 0;
X_rec_new = reshape(X_0, Nx, Ny);

aaa = X_rec_new;
sigma_BP = -1 .* imag((aaa + 1) * eps_b) * w;
eps_BP = real((aaa + 1) * eps_b) ./ const.eps_o;

sigma_BP(zero_index) = 0;
eps_BP(zero_index) = 0;
