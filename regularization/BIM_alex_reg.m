function [dps, ana] = BIM_alex_reg(freq, doi, exp_target_filename, exp_Cal1_filename, exp_Cal2_filename)

pm = path_manager();
pm.add_tomlib();
reg = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/regularizers/tik_matrix_100d_4mm_800MHz.mat");
material_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/system/Material Properties");
cm = DP(freq, material_path + "/cm.txt");
% cal1 = DP(freq, material_path + "/Cal1_Prop.txt");
% cal2 = DP(freq, material_path + "/Cal2_Prop.txt");


%%% -------------------------- Find the Probes_Rx mesh index -------------------------- %%%

tmp_x = ones(doi.total_n, 1);
tmp_y = ones(doi.total_n, 1);

Probes_Rx_coord = zeros(doi.source_num, 2);
Probes_Rx_index = zeros(doi.source_num, 2);

for kk = 1 : source_num
    dist = sqrt((doi.axis_x - Probes_Tx_x(kk) .* tmp_x) .^ 2 + (doi.axis_y - Probes_Tx_y(kk) .* tmp_y) .^ 2);
    [~, min_pos] = min(dist);
    
    [index_x, index_y] = ind2sub([doi.Nx, doi.Ny], min_pos);
    
    Probes_Rx_coord(kk, 1) = doi.axis_x(min_pos);
    Probes_Rx_coord(kk, 2) = doi.axis_y(min_pos);
    
    Probes_Rx_index(kk, 1) = index_x;
    Probes_Rx_index(kk, 2) = index_y;
end
    
%%% ---------------------------- Define the incident field --------------------------- %%%

Ez_inc_MoM = zeros(doi.total_n, doi.src_Tx_N);
Ez_inc_Rx = zeros(doi.src_Rx_N, doi.src_Rx_N);

for src_dash = 1 : doi.src_Tx_N
 
    Probes_pho = sqrt((doi.axis_x(:) - Probes_Tx(src_dash, 1)) .^ 2 + (doi.axis_y(:) - Probes_Tx(src_dash, 2)) .^ 2);
    Ez_inc_MoM(:, src_dash) = -1 * cm.w * cm.uo * 0.25 * besselh(0, 2, cm.kb .* Probes_pho');

    Ez_inc_MoM_rec = reshape(Ez_inc_MoM(:, src_dash), doi.Nx, doi.Ny);

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

eps_tmp(doi.maskp) = eps_r_b;
sigma_tmp(doi.maskp) = sigma_b;
eps_tmp = eps_tmp * eps_o;

X_tmp = (eps_tmp ./ eps_b) - 1i * (sigma_tmp ./ (w * eps_b)) - 1;

%%% ---------------------------- Build the Green function ----------------------------- %%%

al = 4e-3;
a = al / (sqrt(pi));

total_n = length(doi.axis_x);

Gezz = zeros(doi.total_n, doi.total_n);
kb = cm.kb;
for m = 1 : total_n
    x = doi.axis_x(m);
    y = doi.axis_y(m);

    p = sqrt((x - dio.axis_x(:)).^2 + (y - doi.axis_y(:)).^2);
    
     %%% Green function for Ez field  
    Gezz(m, :) = (kb ^ 2) .* ((1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);
    Gezz(m, m) = (kb ^ 2) * ((1 / (kb ^ 2)) + (1i * pi * a * besselh(1, 2, kb * a)) / (2 * kb));
end

Gezz_source = zeros(src_Rx_N, total_n);
for m = 1 : doi.src_Rx_N
    x = Probes_Rx(m, 1);
    y = Probes_Rx(m, 2);

    p = sqrt((x - doi.axis_x).^2 + (y - doi.axis_y).^2);       
    
    %%% Green function for Ez field  
    Gezz_source(m, :) = (kb ^ 2) .* ((-1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);   
end

Gezz = -1 .* Gezz;

%%% ------------------------------ Calibration the measured S-parameters ------------------------------- %%%

Exp_Cal1 = sparameters(exp_Cal1_filename).rfinterp1(freq).Parameters;
Exp_Cal2 = sparameters(exp_Cal2_filename).rfinterp1(freq).Parameters;

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
XX(doi.maskp) = 0;
% XX = X_tmp;
Ez_tot_BIM = zeros(total_n, src_Tx_N);
Lambda = 5e2;

error_obj = zeros(20, 1);
error_tmp = zeros(20, 1);

iter = 20;
XX_array = zeros(doi.Nx, doi.Ny, 2, iter);

tic
for nn = 1 : iter
    
    B = [];
    I = eye(doi.total_n);
    
    GG_Ez_BIM = I + Gezz .* repmat(XX.', doi.total_n, 1);
    
    inv_GG_Ez_BIM = GG_Ez_BIM \ I;
    
    for src_num = 1 : doi.src_Tx_N
        Ez_tot_BIM(:, doi.src_num) = inv_GG_Ez_BIM * Ez_inc_MoM(:, doi.src_num);

        %         fprintf('Finished the %ith source...\n', src_num);
    end
    
    for src_num = 1 : src_Rx_N
        B_dash = Gezz_source .* repmat(Ez_tot_BIM(:, doi.src_num).', doi.src_Rx_N, 1);
        B = [B ; B_dash];
    end
    
    Breg = [B ; Lambda * eye(doi.total_n)];
    
    x0 = Lambda .* X_tmp;
%     x0 = Lambda * zeros(total_n, 1);
    Ereg = [Ess(:) ; x0];
    
%     XX = lsqr(Breg, Ereg, 1e-6, 5000);
%     XX = lsmr(B, Ess(:), Lambda);
    XX = lsmr(Breg, Ereg);
    
    XX(doi.maskp) = 0;
    
%     error_obj(nn) = norm(B * XX - Ess(:)) ^ 2 / norm(Ess(:)) ^ 2;
    error_obj(nn) = sqrt(sum(abs(B * XX - Ess(:)) .^ 2) ./ sum(abs(Ess(:)) .^ 2));
%     error_tmp(nn) = norm(XX - X_tmp) ^ 2 / norm(X_tmp) ^ 2;
    XX_rmse = XX;
    XX_rmse(doi.maskp) = 1;
    XX_tmp_rmse = X_tmp;
    XX_tmp_rmse(doi.maskp) = 1;
    error_tmp(nn) = sqrt(sum(abs(XX_rmse - XX_tmp_rmse) .^ 2) ./ sum(abs(XX_tmp_rmse) .^ 2));
    
    fprintf('The %ith iteration\n', nn);
    fprintf('...Error: %i\n', error_obj(nn));

    XX_rec = reshape(XX, doi.Nx, doi.Ny);
    complex_eps = (XX_rec + 1) * cm.eps_complex;
    eps_r = real(complex_eps ./ cm.eps_o);
    sig = -1 .* imag(complex_eps * cm.w);
    XX_array(:, :, 1, nn) = eps_r;
    XX_array(:, :, 2, nn) = sig;
    

end
toc

dps = DD(eps_r, sig, freq);
ana.error_obj = error_obj;
ana.error_tmp = error_tmp;
ana.XX_array = XX_array;

fprintf('The %ith frequency sample is finished...\n', freq_num)    
end