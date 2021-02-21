function green = MoMGreenFunc(stp, Ez_inc_MoM)

% ======================= Attempting to Load a saved version ==============
pm = path_manager();
green_path = pm.data + "/cache/MoMGreen_" + string(stp.freq) + "Hz.mat";
green_path = char(green_path);  % path could be a symlink from home to somewhere else.
try
    green = load(green_path);
    green = green.green;
    fprintf('Green function loaded from cache...\n');
    fprintf(green_path)
    return  % skip the reaminder of the function.
catch
    % continue executing the code body below:
end


%%% Build the model for calibration phantoms (Contrasts)----------------------------- %%%

X_Cal1 = buildCalModel(stp, stp.cal1);
X_Cal2 = buildCalModel(stp, stp.cal2);

%%% ---------------------------- Build the Green function ----------------------------- %%%

al = stp.res * 1e-3;
a = al / (sqrt(pi));
Gezz = zeros(stp.total_n, stp.total_n);
kb = stp.bg.kb;
for m = 1 : stp.total_n
    x = stp.axis_x(m);
    y = stp.axis_y(m);

    p = sqrt((x - stp.axis_x(:)).^2 + (y - stp.axis_y(:)).^2);

     %%% Green function for Ez field  
    Gezz(m, :) = (kb ^ 2) .* ((1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);
    Gezz(m, m) = (kb ^ 2) * ((1 / (kb ^ 2)) + (1i * pi * a * besselh(1, 2, kb * a)) / (2 * kb));
end
Gezz_source = zeros(stp.src_Rx_N, stp.total_n);
for m = 1 : stp.src_Rx_N
    Probes_Rx = stp.probes_Tx;
    x = Probes_Rx(m, 1);
    y = Probes_Rx(m, 2);

    p = sqrt((x - stp.axis_x(:)).^2 + (y - stp.axis_y(:)).^2);       
    
    % singular_index = p <= 1e-5;
    %%% Green function for Ez field  
    Gezz_source(m, :) = (kb ^ 2) .* ((-1i * pi .* a) / (2 * kb)) .* besselj(1, kb .* a) .* besselh(0, 2, kb .* p);   
    Gezz_source(m, p <= 1e-5) = (kb ^ 2) * ((1 / (kb ^ 2)) + (1i * pi * a * besselh(1, 2, kb * a)) / (2 * kb));
end
clear kb
I = eye(stp.total_n);

GG_Ez_Cal1 = I + Gezz .* repmat(X_Cal1.', stp.total_n, 1);
Ez_tot_Cal1 = zeros(stp.total_n, stp.src_Tx_N);

GG_Ez_Cal2 = I + Gezz .* repmat(X_Cal2.', stp.total_n, 1);
Ez_tot_Cal2 = zeros(stp.total_n, stp.src_Tx_N);

for src_num = 1 : stp.src_Tx_N

    Ez_tot_Cal1(:, src_num) = GG_Ez_Cal1 \ Ez_inc_MoM(:, src_num);
    Ez_tot_Cal2(:, src_num) = GG_Ez_Cal2 \ Ez_inc_MoM(:, src_num);

end

wr_Ez_Cal1 = repmat(X_Cal1, 1, stp.src_Tx_N) .* Ez_tot_Cal1;
Ess_Cal1 = Gezz_source * wr_Ez_Cal1;
wr_Ez_Cal2 = repmat(X_Cal2, 1, stp.src_Tx_N) .* Ez_tot_Cal2;
Ess_Cal2 = Gezz_source * wr_Ez_Cal2;

green.Ess_Cal1 = Ess_Cal1;
green.Ess_Cal2 = Ess_Cal2;
green.Gezz = -1 * Gezz;
green.Gezz_source = Gezz_source;

if exist(pm.data + "/cache", "dir") ~= 0
   save(green_path, "green");
   fprintf("Green functions saved @ " + string(green_path));
end
fprintf('Calculation of Green functions and MoM data for Cal1 and Cal2 are finished...\n');
end

