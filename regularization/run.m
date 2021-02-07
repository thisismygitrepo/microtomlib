
clc
clear

pm = path_manager();

% database_path = fullfile(home, "emvision/Algorithm/toml_data/solvers/Parameters_for_solvers");
file_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/large_stroke_simulations/example1/wIXI087-Guys-0768-T1_t0_1.s16p");
cal1_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/cals/Cal1.s16p");
cal2_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/cals/Cal2.s16p");


% before running this, make sure your current directory is @ toml/solvers
[eps, sig] = BIM_alex_reg(700e6, file_path, cal1_path, cal2_path);

