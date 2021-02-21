
% This code is only used for testing and debug the BP solver

clc
clear

pm = path_manager();
sources = ["simulation", "PA"];  % Which data?
source = sources(1);  % <---------- CHANGE THIS TO (1 OR 2).

if source == "simulation"
    target_filename = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/simulations\wiener\large_stroke_simulations\example1\wIXI087-Guys-0768-T1_t0_1.s16p");
    cal1_filename = pm.join(pm.gdrive, "emvision/Algorithm/toml_data\simulations\wiener\cals\Cal1.s16p");
    cal2_filename = pm.join(pm.gdrive, "emvision/Algorithm/toml_data\simulations\wiener\cals\Cal2.s16p");
    incident_filename = pm.join(pm.gdrive, "emvision\Algorithm\toml_data\simulations\wiener\cals\S_Incident.s16p");
elseif source == "PA"
    incident_filename = pm.join(pm.gdrive, "emvision\Algorithm\toml_data\simulations\wiener\cals\S_Incident.s16p");    
    target_filename = pm.join(pm.gdrive, "emvdata/MeasurementData/Clinical-Data-20200220-Stage-1/output/000000/2-delta+0d0h8m/dscan1-emscan-e00-cal.s16p");
    cal1_filename = pm.join(pm.gdrive, "emvdata/MeasurementData/Clinical-Data-20200220-Stage-1/output/000000/2-delta+0d0h8m/cal/dscan1-cal1-cal.s16p");
    cal2_filename = pm.join(pm.gdrive, "emvdata/MeasurementData/Clinical-Data-20200220-Stage-1/output/000000/2-delta+0d0h8m/cal/dscan1-cal2-cal.s16p");
end


resolution = 4;  % in mm.
freq = 700e6;  % in Hz.
pm.add_tomlib()
material_path = pm.join(pm.gdrive, "emvision/Algorithm/toml_data/system/Material Properties");
cm = DP(freq, material_path + "/cm.txt");
cal1 = DP(freq, material_path + "/Cal1_Prop.txt");
cal2 = DP(freq, material_path + "/Cal2_Prop.txt");
stp = SystemSetup(resolution, cm, cal1, cal2);

tic
props = BP(stp, target_filename, cal1_filename, cal2_filename, incident_filename); 
mkdir results
save("results/" + source + "_" + "results.mat", "props")
toc


% Compare Sim to PA.

% pa = load("results/simulation_results.mat");
% sim = load("results/PA_results.mat");
% bias = (sim.eps_BP == 0) * 50;
% im_pa = imagesc(pa.eps_BP + bias);
% colorbar
% figure
% im_sim = imagesc(sim.eps_BP + bias);
% colorbar

