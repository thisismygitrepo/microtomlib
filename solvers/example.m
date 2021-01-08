clc
clear

pm = path_manager();
database_path = fullfile(pm.gdrive, "emvision/Algorithm/toml_data/solvers/Parameters_for_solvers");
file_path = fullfile(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/large_stroke_simulations/example1/wIXI087-Guys-0768-T1_t0_1.s16p");
cal1_path = fullfile(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/cals/Cal1.s16p");
cal2_path = fullfile(pm.gdrive, "emvision/Algorithm/toml_data/simulations/wiener/cals/Cal2.s16p");


% before running this, make sure your current directory is @ toml/solvers
[eps, sig] = run_multiFreq_solvers("MRCSI", file_path, cal1_path, cal2_path, database_path);
