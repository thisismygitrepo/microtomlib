function dps = BP(stp, target_filename, Cal1_filename, Cal2_filename, Incident_filename)


%%% ---------------------------- Define the incident field --------------------------- %%%

bg = stp.bg;
Ez_inc_MoM = zeros(stp.total_n, stp.src_Tx_N);
for src_dash = 1 : stp.src_Tx_N
    Probes_pho = sqrt((stp.axis_x(:) - stp.probes_Tx(src_dash, 1)) .^ 2 + (stp.axis_y(:) - stp.probes_Tx(src_dash, 2)) .^ 2);
    Ez_inc_MoM(:, src_dash) = -1 * bg.w * bg.uo * 0.25 * besselh(0, 2, bg.kb .* Probes_pho');
end
green = MoMGreenFunc(stp, Ez_inc_MoM);


%%% -- Do the calibration (Match 3D SParams to 2D Field) ------------------ %%%

CST_Cal1 = sparameters(Cal1_filename).rfinterp1(stp.freq).Parameters;
CST_Cal2 = sparameters(Cal2_filename).rfinterp1(stp.freq).Parameters;
CST_case = sparameters(target_filename).rfinterp1(stp.freq).Parameters;
if nargin == 6
    CST_Inc = sparameters(Incident_filename).rfinterp1(stp.freq).Parameters;
    CST_Cal1 = CST_Cal1 - CST_Inc;
    CST_Cal2 = CST_Cal2 - CST_Inc;
    CST_case = CST_case - CST_Inc;
end
alpha_Ez = (green.Ess_Cal1 - green.Ess_Cal2) ./ (CST_Cal1 - CST_Cal2);
belta_Ez = green.Ess_Cal1 - alpha_Ez .* CST_Cal1;
S_scat = alpha_Ez .* CST_case + belta_Ez;


%%% ---------------------------- The BP solver --------------------------- %%%

wr_bp_0 = zeros(stp.total_n, stp.src_Tx_N);
u_0 = zeros(stp.total_n, stp.src_Tx_N);
Ge_star_fi = zeros(stp.total_n, stp.src_Tx_N);
Ge_Ge_star_fi = zeros(stp.src_Rx_N, stp.src_Tx_N);

for src = 1 : stp.src_Tx_N
    Ge_star_fi(:, src) = green.Gezz_source' * S_scat(:, src);
    Ge_Ge_star_fi(:, src) = green.Gezz_source * Ge_star_fi(:, src);
    wr_bp_0(:, src) = (norm(Ge_star_fi(:, src)) ^ 2 / norm(Ge_Ge_star_fi(:, src)) ^ 2) * Ge_star_fi(:, src);
    u_0(:, src) = Ez_inc_MoM(:, src) + green.Gezz * wr_bp_0(:, src);
end

X_0 = sum(wr_bp_0 .* conj(u_0), 2) ./ sum(abs(u_0) .^ 2, 2);
X_rec_new = reshape(X_0, stp.Nx, stp.Ny);
X_rec_new(stp.maskp) = 0;
eps_complex = (X_rec_new + 1) * bg.eps_complex;
eps_r = real(eps_complex) ./ bg.eps_o;
sig = -1 .* imag(eps_complex) * bg.w;

pm = path_manager();
pm.add_tomlib();  % We need DD class for output.
dps = DD(eps_r, sig, stp);

end
