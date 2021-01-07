%%% This code is only used for testing and debug the new BP solver

clc
clear

pm = path_manager();

target_filename = pm.gdrive + "\emvision\Algorithm\toml_data\simulations\wiener\large_stroke_simulations\example1\wIXI087-Guys-0768-T1_t0_1.s16p";
Cal1_filename = pm.gdrive + "\emvision\Algorithm\toml_data\simulations\wiener\cals\Cal1.s16p";
Cal2_filename = pm.gdrive + "\emvision\Algorithm\toml_data\simulations\wiener\cals\Cal2.s16p";

[eps_BP, sigma_BP] = func_BP_CST_dynamic(730.4e6, target_filename, Cal1_filename, Cal2_filename, 2);    