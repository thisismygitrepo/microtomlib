function green = MoM_GreenFunc(freq, DOI_setting, Ez_inc_MoM, eps_r_b, sigma_b, Cal1_eps, Cal1_sigma, Cal2_eps, Cal2_sigma, source_num)

eps_o = 8.854187817e-12;
uo = 4e-7 * pi;  

%%% --------------------------- Define the position of the waveguide ----------------------------%%%

Probes_Tx_x = DOI_setting.ant_xy(:, 2);
Probes_Tx_y = DOI_setting.ant_xy(:, 1);

Probes_Tx_x = Probes_Tx_x * 1e-3;
Probes_Tx_y = Probes_Tx_y * 1e-3;

Probes_Tx = zeros(source_num, 2);
Probes_Tx(:, 1) = Probes_Tx_x;
Probes_Tx(:, 2) = Probes_Tx_y;

Probes_Rx = Probes_Tx;

%%% --------------------------- Define the computation domain ----------------------------- %%%

src_Tx_N = length(Probes_Tx);
src_Rx_N = length(Probes_Rx);

x_dash = DOI_setting.y_axis;
y_dash = DOI_setting.x_axis;

x_dash = x_dash * 1e-3;
y_dash = y_dash * 1e-3;

Nx = length(x_dash);
Ny = length(y_dash);

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

Axis_x = Axis_x(:);
Axis_y = Axis_y(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 2 * pi * freq;
eps_b = eps_r_b * eps_o - 1i * sigma_b / w;
kb = w * sqrt(uo * eps_b);

%%% ---------------------------- Build the model for calibration phantoms ----------------------------- %%%

%%% Build the contrast for Cal1
[eps_gama_Cal1, sigma_gama_Cal1] = build_Exp_Cal_model(x_dash, y_dash, Cal1_eps, Cal1_sigma, eps_r_b, sigma_b);

eps_gama_Cal1 = eps_gama_Cal1 .* eps_o;
X_Cal1 = (eps_gama_Cal1 ./ eps_b) - 1i * (sigma_gama_Cal1 ./ (w * eps_b)) - 1; 
X_Cal1 = X_Cal1(:);

%%% Build the contrast for Cal2
[eps_gama_Cal2, sigma_gama_Cal2] = build_Exp_Cal_model(x_dash, y_dash, Cal2_eps, Cal2_sigma, eps_r_b, sigma_b);

eps_gama_Cal2 = eps_gama_Cal2 .* eps_o;
X_Cal2 = (eps_gama_Cal2 ./ eps_b) - 1i * (sigma_gama_Cal2 ./ (w * eps_b)) - 1; 
X_Cal2 = X_Cal2(:);

%%% ---------------------------- Build the Green function ----------------------------- %%%

al = DOI_setting.res * 1e-3;
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
    
    singular_index = find(p <= 1e-5);
    %%% Green function for Ez field  
    Gezz_source(m, :) = (kb ^ 2) .* ((-1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);   
    Gezz_source(m, singular_index) = (kb ^ 2) * ((1 / (kb ^ 2)) + (1i * pi * a * besselh(1, 2, kb * a)) / (2 * kb));
end

I = eye(total_n);

GG_Ez_Cal1 = I + Gezz .* repmat(X_Cal1.', total_n, 1);
Ez_tot_Cal1 = zeros(total_n, src_Tx_N);

GG_Ez_Cal2 = I + Gezz .* repmat(X_Cal2.', total_n, 1);
Ez_tot_Cal2 = zeros(total_n, src_Tx_N);

for src_num = 1 : src_Tx_N

    Ez_tot_Cal1(:, src_num) = GG_Ez_Cal1 \ Ez_inc_MoM(:, src_num);
    Ez_tot_Cal2(:, src_num) = GG_Ez_Cal2 \ Ez_inc_MoM(:, src_num);

end

wr_Ez_Cal1 = repmat(X_Cal1, 1, src_Tx_N) .* Ez_tot_Cal1;
Ess_Cal1 = Gezz_source * wr_Ez_Cal1;

wr_Ez_Cal2 = repmat(X_Cal2, 1, src_Tx_N) .* Ez_tot_Cal2;
Ess_Cal2 = Gezz_source * wr_Ez_Cal2;

green.Ess_Cal1 = Ess_Cal1;
green.Ess_Cal2 = Ess_Cal2;
green.Gezz = Gezz;
green.Gezz_source = Gezz_source;

fprintf('Calculation of Green functions and MoM data for Cal1 and Cal2 are finished...\n');


