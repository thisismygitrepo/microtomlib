function [eps_BP, sigma_BP] = func_BP_CST(freq_num, target_filename, resolution)

freq = [701; 750.5; 800; 851; 900.5; 950; 1001; 1050.5; 1100; 1151; 1200.5; 1250; 1301] .* 1e6;
w = 2 * pi * freq(freq_num);
eps_o = 8.854187817e-12;
uo = 4e-7 * pi;  
c = 299792458;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The corresponding permittivity in eps_r_b is for 700MHz, 750MHz, 800MHz, and 850MHz, and the same for sigma_b

eps_r_b = [52.4; 52.06; 52.665; 52.203; 52.1464; 52.1354; 52.8048; 52.6253; 52.0931; 52.4981; 52.102; 51.9399; 51.7879];
sigma_b = [2.9928; 2.5472; 2.7799; 3.3259; 2.9245; 3.2465; 3.0125; 3.9947; 3.859; 3.9225; 4.3962; 4.1691; 4.6537] .* w .* eps_o;
eps_b = eps_r_b(freq_num) * eps_o - 1i * sigma_b(freq_num) / w;
kb = w * sqrt(uo * eps_b);
eps_bx = eps_r_b(freq_num) * eps_o - 1i * sigma_b(freq_num) / w;

source_num = 16;

%%% --------------------------- Define the position of the waveguide ----------------------------%%%

Probes_Tx_x = [106.4; 83.4; 51; 17.23; -16.94; -51.73; -84.34; -106.2; -106.2; -84.34; -51.73; -16.94; 17.23; 51; 83.4; 106.4] .* 1e-3;
Probes_Tx_y = [21.12; 55.56; 75.96; 85.36; 86.43; 77.79; 56.53; 21.18; -21.18; -56.53; -77.8; -86.43; -85.36; -75.96; -55.56; -21.12] .* 1e-3;

Probes_Tx = zeros(source_num, 2);
Probes_Tx(:, 1) = Probes_Tx_x;
Probes_Tx(:, 2) = Probes_Tx_y;

Probes_Rx = Probes_Tx;

src_Tx_N = length(Probes_Tx(:, 1));
src_Rx_N = length(Probes_Rx(:, 1));

%%% --------------------------- Define the computation domain ----------------------------- %%%

x_dash = [-110 : 4 : 110] .* 1e-3;
y_dash = [-93 : 4 : 92] .* 1e-3;

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

for src_dash = 1 : src_Tx_N
 
    Probes_pho = sqrt((Axis_x(:) - Probes_Tx(src_dash, 1)) .^ 2 + (Axis_y(:) - Probes_Tx(src_dash, 2)) .^ 2);
    Ez_inc_MoM(:, src_dash) = -1 * w * uo * 0.25 * besselh(0, 2, kb .* Probes_pho');
    
end

%%% ---------------------------- Build the Green function ----------------------------- %%%

al = 4e-3;
a = al / (sqrt(pi));

total_n = length(Axis_x);


I = eye(total_n);

if freq_num == 1
    load Gezz_700MHz Gezz
    load Gezz_source_700MHz Gezz_source
elseif freq_num == 2
    load Gezz_750MHz Gezz
    load Gezz_source_750MHz Gezz_source
elseif freq_num == 3
    load Gezz_800MHz Gezz
    load Gezz_source_800MHz Gezz_source
elseif freq_num == 4
    load Gezz_850MHz Gezz
    load Gezz_source_850MHz Gezz_source
elseif freq_num == 5
    load Gezz_900MHz Gezz
    load Gezz_source_900MHz Gezz_source
elseif freq_num == 6
    load Gezz_950MHz Gezz
    load Gezz_source_950MHz Gezz_source
elseif freq_num == 7
    load Gezz_1000MHz Gezz
    load Gezz_source_1000MHz Gezz_source
elseif freq_num == 8
    load Gezz_1050MHz Gezz
    load Gezz_source_1050MHz Gezz_source
elseif freq_num == 9
    load Gezz_1100MHz Gezz
    load Gezz_source_1100MHz Gezz_source
elseif freq_num == 10
    load Gezz_1150MHz Gezz
    load Gezz_source_1150MHz Gezz_source
elseif freq_num == 11
    load Gezz_1200MHz Gezz
    load Gezz_source_1200MHz Gezz_source
elseif freq_num == 12
    load Gezz_1250MHz Gezz
    load Gezz_source_1250MHz Gezz_source
elseif freq_num == 13
    load Gezz_1300MHz Gezz
    load Gezz_source_1300MHz Gezz_source
end

%%% ------------------------------ Calibration the measured S-parameters ------------------------------- %%%

% addpath('C:\MATLAB\HeadScanner\tomography-integral\Machine Learning\Alex');

CST_Cal1_struct = sparameters('Cal1.s16p');
CST_Cal2_struct = sparameters('Cal2.s16p');
CST_inc_struct = sparameters('S_Incident.s16p');

CST_Cal1 = -1 .* CST_Cal1_struct.Parameters;
CST_Cal2 = -1 .* CST_Cal2_struct.Parameters;
CST_inc = -1 .* CST_inc_struct.Parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% The freq_sample corresponds to 700MHz, 750MHz, 800MHz,..., 1300MHz 
freq_sample = [135; 168; 201; 235; 268; 301; 335; 368; 401; 435; 468; 501; 535];

CST_Cal1 = CST_Cal1(:, :, freq_sample(freq_num));
CST_Cal2 = CST_Cal2(:, :, freq_sample(freq_num));
CST_inc = CST_inc(:, :, freq_sample(freq_num));

CST_Cal1_scat = CST_Cal1 - CST_inc;
CST_Cal2_scat = CST_Cal2 - CST_inc;

load Ess_MoM_Cal1_array.mat Ess_MoM_Cal1_array
load Ess_MoM_Cal2_array.mat Ess_MoM_Cal2_array
Ess_MoM_Cal1 = Ess_MoM_Cal1_array(:, :, freq_num);
Ess_MoM_Cal2 = Ess_MoM_Cal2_array(:, :, freq_num);

alpha_Ez = (Ess_MoM_Cal1 - Ess_MoM_Cal2) ./ (CST_Cal1_scat - CST_Cal2_scat);
belta_Ez = Ess_MoM_Cal1 - alpha_Ez .* CST_Cal1_scat;

CST_case_struct = sparameters(target_filename);
CST_case = -1 .* CST_case_struct.Parameters;
CST_case = CST_case(:, :, freq_sample(freq_num));
CST_case = CST_case - CST_inc;
% CST_case = CST_inc - CST_inc;

S_scat = alpha_Ez .* CST_case + belta_Ez;
% Ess = Ess - Ez_inc_Rx;

%%% ----------------------------------------- The CSI inverse solver ------------------------------------------- %%%

load CST_model_mask_4mm.mat CST_model_mask
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

switch resolution 
    case '4mm'
        aaa = X_rec_new;
        sigma_BP = -1 .* imag((aaa + 1) * eps_bx) * w;
        eps_BP = real((aaa + 1) * eps_bx) ./ eps_o;
        
        sigma_BP(zero_index) = 0;
        eps_BP(zero_index) = 0;
    case '2mm'
        aaa = X_rec_new;
        sigma_BP = -1 .* imag((aaa + 1) * eps_bx) * w;
        eps_BP = real((aaa + 1) * eps_bx) ./ eps_o;
        
        sigma_BP(zero_index) = 0;
        eps_BP(zero_index) = 0;
        
        sigma_BP = imresize(sigma_BP, 2);
        eps_BP = imresize(eps_BP, 2);
        
        eps_BP(eps_BP < 0) = 0;
        sigma_BP(sigma_BP < 0) = 0;
end
% sigma_BP(zero_index) = 0;
% eps_BP(zero_index) = 0;
% figure; imagesc(eps_csi);axis image;axis off;colormap(jet);caxis([35 60])
% figure; imagesc(sigma_csi);axis image;axis off;colormap(jet);caxis([0 0.3])
% pause(0.001)


















