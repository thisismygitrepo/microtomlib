function [eps_csi, sigma_csi] = MR_CSI_newReg_v1(freq_num, exp_target_filename, exp_Cal1_filename, exp_Cal2_filename, database_path, itr_num)

%%% The DOI has the physical range of x: -100mm : 100mm, y: -115mm : 115mm
load(database_path + "/doi_masks.mat")

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

dx = 4e-3;
dy = 4e-3;

x_dash = [-115e-3 : dx : 115e-3];
y_dash = [-100e-3 : dy : 96e-3];

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

freq_sample = [135 : 4 : 467];     %%% Change this to [135 : 4 : 467] for Alex's CST simulations

Exp_Cal1 = Exp_Cal1(:, :, freq_sample(freq_num));
Exp_Cal2 = Exp_Cal2(:, :, freq_sample(freq_num));

Et_MoM_Cal1 = Ez_tot_MoM_array_Cal1(:, :, freq_num);
Et_MoM_Cal2 = Ez_tot_MoM_array_Cal2(:, :, freq_num);

Ess_MoM_Cal1 = Et_MoM_Cal1 - Ez_inc_Rx;
Ess_MoM_Cal2 = Et_MoM_Cal2 - Ez_inc_Rx;

% alpha_Ez = (Ess_MoM_Cal1 - Ess_MoM_Cal2) ./ (Exp_Cal1 - Exp_Cal2);
% belta_Ez = Ess_MoM_Cal1 - alpha_Ez .* Exp_Cal1;

% Exp_case = importdata(exp_target_filename);
% Exp_case = Exp_case(:, :, freq_num);

Exp_case_struct = sparameters(exp_target_filename);
Exp_case = Exp_case_struct.Parameters;
Exp_case = Exp_case(:, :, freq_sample(freq_num));

% Exp_case = Exp_case - Exp_Cal1;

% S_scat = alpha_Ez .* Exp_case + belta_Ez;

alpha_Ez = Ess_MoM_Cal2 ./ Exp_Cal2;
S_scat = alpha_Ez .* Exp_case;


%%% ----------------------------------------------------- The CSI inverse solver -----------------------------------------------------

max_itr = 10000;
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

X_0(stp.maskp) = 0;

for kk = 1 : src_Tx_N
    wr_bp_0(stp.maskp, kk) = 0;
end

%%% ---------------------------------- CSI (upgrade the contrast source)--------------------------------
eta_d = zeros(max_itr, 1);
g_wr = zeros(total_n, src_Tx_N);                   % The valua of g_wr in the first iteration is not used, so it is zero
g_x = zeros(total_n);                           % The value of g_x in the first iteration is not used, so it is zero  
v = zeros(total_n, src_Tx_N);                   % The value of v in the first iteration is set to be zero
d = zeros(total_n);                          % The value of d in the first iteration is set to be zero
alpha_wr = zeros(max_itr, 1);                            % The first value of alpha_wr is not used so it is zero
alpha_x = zeros(max_itr, 1);                             % The first value of alpha_x is not used so it is zero
gama_wr = zeros(max_itr, 1);                                % The values of gama_wr in the first and second iteration are not used, so they are zeros
gama_x = zeros(max_itr, 1);
Fd = zeros(max_itr, 1);                                     
Fs = zeros(max_itr, 1);                                     % The first value of Fs will not be used in the algorithm, so we can set it as zero
F_TV = zeros(max_itr, 1);


%%% Initialization
X_csi = X_0;
wr_Ez_csi = wr_bp_0;
pho = S_scat - (Gezz_source * wr_bp_0);
Ez_tot_csi = Ez_inc_MoM + Gezz * wr_bp_0;
r = repmat(X_csi, 1, src_Tx_N) .* Ez_tot_csi - wr_bp_0;
eta_s = 1 / (square_norm_sum(S_scat, dx, dy));
eta_d(1) = 1 / (square_norm_sum(repmat(X_csi, 1, src_Tx_N) .* Ez_inc_MoM, dx, dy));
Fd(1) = eta_d(1) * square_norm_sum((repmat(X_csi, 1, src_Tx_N) .* Ez_tot_csi) - wr_Ez_csi, dx, dy);

X_rec = reshape(X_csi, Nx, Ny);
[X_gy, X_gx] = gradient(X_rec, y_dash, x_dash);
% [X_gy, X_gx] = gradient(X_rec);
abs_gradientX_square = X_gx(:) .* conj(X_gx(:)) + X_gy(:) .* conj(X_gy(:));

for itr = 2 : itr_num
    %%% Update the contrast source (wr) and the total E field
    Ge_star_pho = (Gezz_source' * pho);
    
    X_csi_adj = conj(X_csi);
    Ge_star_X_r = Gezz' * (repmat(X_csi_adj, 1, src_Tx_N) .* r);
    
    g_wr_dash = g_wr;
    
    g_wr = -1 * eta_s .* Ge_star_pho - eta_d(itr - 1) .* (r - Ge_star_X_r);
    
    if itr == 2
        v = g_wr;
    else
        numerator = real(sum(inner_product(g_wr, g_wr - g_wr_dash, dx, dy)));
        denominator = sum(inner_product(g_wr_dash, g_wr_dash, dx, dy));
        gama_wr(itr) = numerator / denominator;
        v = g_wr + gama_wr(itr) .* v;
    end
    
    Ge_source_v= Gezz_source * v;
    Ge_v = Gezz * v;
    
    numerator = -1 * real(sum(inner_product(g_wr, v, dx, dy)));
    denominator = eta_s * square_norm_sum(Ge_source_v, dx, dy) + eta_d(itr - 1) * square_norm_sum(v - repmat(X_csi, 1, src_Tx_N) .* Ge_v, dx, dy);
    alpha_wr(itr) = numerator / denominator;
    
    wr_Ez_csi = wr_Ez_csi + alpha_wr(itr) .* v;
    Ez_tot_csi = Ez_inc_MoM + Gezz * wr_Ez_csi;
    
    for kk = 1 : src_Tx_N
        wr_Ez_csi(stp.maskp, kk) = 0;
    end
    
    %%% Update the contrast (X)
    Fd(itr) = square_norm_sum((repmat(X_csi, 1, src_Tx_N) .* Ez_tot_csi) - wr_Ez_csi, dx, dy);
    Fs(itr) = square_norm_sum(S_scat - Gezz_source * wr_Ez_csi, dx, dy) ./ square_norm_sum(S_scat, dx, dy);
    
    %%% ------------------------------------- Add the multiplicative TV reguarization ----------------------------------------- %%%
    mesh_size = dx;
    delta_square = Fd(itr) * (1 / mesh_size) ^ 2;
    
    V = (dx ^ 2) * Nx * Ny;
    
    X_rec = reshape(X_csi, Nx, Ny);
    [X_gy, X_gx] = gradient(X_rec, y_dash, x_dash);
%     [X_gy, X_gx] = gradient(X_rec);
    
    abs_gradientX_square_dash = abs_gradientX_square;
    abs_gradientX_square = X_gx(:) .* conj(X_gx(:)) + X_gy(:) .* conj(X_gy(:));
    
    g_TV_tmp_x = X_gx(:) ./ (abs_gradientX_square + delta_square);
    g_TV_tmp_x = reshape(g_TV_tmp_x, Nx, Ny);
    g_TV_tmp_y = X_gy(:) ./ (abs_gradientX_square + delta_square);
    g_TV_tmp_y = reshape(g_TV_tmp_y, Nx, Ny);
    
    g_TV_rec = (1 / V) .* divergence(y_dash, x_dash, g_TV_tmp_y, g_TV_tmp_x);
%     g_TV_rec = (1 / V) .* divergence(g_TV_tmp_y, g_TV_tmp_x);
    g_TV_rec = g_TV_rec(:);
    g_TV = g_TV_rec;
    
    numerator_left = eta_d(itr - 1) .* sum((wr_Ez_csi - repmat(X_csi, 1, src_Tx_N) .* Ez_tot_csi) .* conj(Ez_tot_csi), 2);
    numerator_right = (Fd(itr) + Fs(itr)) .* g_TV;
    denominator = sum(abs(Ez_tot_csi) .^ 2, 2);
    
    g_x_dash = g_x;
    
    g_x = (numerator_left + numerator_right) ./ denominator;
    
    if itr == 2
       d = g_x;
    else
        numerator = real(inner_product(g_x, g_x - g_x_dash, dx, dy));
        denominator = inner_product(g_x_dash, g_x_dash, dx, dy);
        gama_x(itr) = numerator / denominator;
        d = g_x + gama_x(itr) .* d;
    end
    
    A_alpha = Fs(itr);
    B_alpha = Fd(itr);
    C_alpha = 2 * eta_d(itr - 1) * real(sum(inner_product(repmat(d, 1, src_Tx_N) .* Ez_tot_csi, repmat(X_csi, 1, src_Tx_N) .* Ez_tot_csi - wr_Ez_csi, dx, dy)));
    D_alpha = eta_d(itr - 1) * square_norm_sum(repmat(d, 1, src_Tx_N) .* Ez_tot_csi, dx, dy);
    
    b_alpha = (1 / sqrt(V)) .* (1 ./ sqrt(abs_gradientX_square + delta_square));
    d_rec = reshape(d, Nx, Ny);
    [d_gy, d_gx] = gradient(d_rec, y_dash, x_dash);
%     [d_gy, d_gx] = gradient(d_rec);

    E_alpha = 2 * real(gradient_inner_product(b_alpha .* X_gx(:), b_alpha .* X_gy(:), b_alpha .* d_gx(:), b_alpha .* d_gy(:), dx, dy));
    F_alpha = gradient_inner_product(b_alpha .* d_gx(:), b_alpha .* d_gy(:), b_alpha .* d_gx(:), b_alpha .* d_gy(:), dx, dy);
    G_alpha = 1;
    
    poly = [4 * D_alpha * F_alpha; 3 * C_alpha * F_alpha + 3 * D_alpha * E_alpha; 2 * C_alpha * E_alpha + 2 * D_alpha * G_alpha + 2 * A_alpha * F_alpha + 2 * B_alpha * F_alpha; C_alpha * G_alpha + A_alpha * E_alpha + B_alpha * E_alpha];
    alpha_root = roots(poly);
    
    for count = 1 : 3
        if isreal(alpha_root(count))
            alpha_x(itr) = alpha_root(count);
        end
    end
    
    %%% ------------------------------------- Update contrast X_csi using MR-TV ------------------------------------- %%%
    
    X_csi = X_csi + alpha_x(itr) .* d; 
    X_csi(stp.maskp) = 0;
    
    %%% -------------------------------------------------------------------------------------------------------------------

    
    %%% Update the parameters in CSI (pho, r, eta_d)
    pho = S_scat - (Gezz_source * wr_Ez_csi);
    r = repmat(X_csi, 1, src_Tx_N) .* Ez_tot_csi - wr_Ez_csi;
    eta_d(itr) = 1 / (square_norm_sum(repmat(X_csi, 1, src_Tx_N) .* Ez_inc_MoM, dx, dy));
    
    F_TV(itr) = sum((abs_gradientX_square + delta_square) ./ (abs_gradientX_square_dash + delta_square) * al ^ 2) ./ V;

    X_rec_new = reshape(X_csi, Nx, Ny);
    
    sigma_csi = -1 .* imag((X_rec_new + 1) * eps_b) * w;
    eps_csi = real((X_rec_new + 1) * eps_b) ./ eps_o;
    
    
end

 fprintf('The %ith frequency sample is finished...\n', freq_num)   
