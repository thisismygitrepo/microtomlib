
%%% This code is only used for testing and debug the new BP solver

clc
clear

sources = ["simulation", "PA"];
source = sources(1);  % <---------- CHANGE THIS TO (1 OR 2).

pm = path_manager();


if source == "simulation"

    target_filename = pm.join(pm.gdrive, "emvision\Platform\System Simulations\Alex\Results\BP_Solver_Lei\105014_1.s16p");
    Cal1_filename = pm.join(pm.gdrive, "emvision\Platform\System Simulations\Alex\Results\BP_Solver_Lei\Parameters\Cal1.s16p");
    Cal2_filename = pm.join(pm.gdrive, "emvision\Platform\System Simulations\Alex\Results\BP_Solver_Lei\Parameters\Cal2.s16p");
    Incident_filename = pm.join(pm.gdrive, "emvision\Algorithm\toml_data\simulations\wiener\cals\S_Incident.s16p");

elseif source == "PA"

    target_filename = pm.join(pm.gdrive, "emvdata/MeasurementData/Clinical-Data-20200220-Stage-1/output/000000/2-delta+0d0h8m/dscan1-emscan-e00-cal.s16p");
    Cal1_filename = pm.join(pm.gdrive, "emvdata/MeasurementData/Clinical-Data-20200220-Stage-1/output/000000/2-delta+0d0h8m/cal/dscan1-cal1-cal.s16p");
    Cal2_filename = pm.join(pm.gdrive, "emvdata/MeasurementData/Clinical-Data-20200220-Stage-1/output/000000/2-delta+0d0h8m/cal/dscan1-cal2-cal.s16p");
    
end


tic

if source == 'simulation'
    [eps_BP, sigma_BP] = func_BP_CST_dynamic(700e6, 4, target_filename, Cal1_filename, Cal2_filename, Incident_filename); 
else
    [eps_BP, sigma_BP] = func_BP_CST_dynamic(700e6, 4, target_filename, Cal1_filename, Cal2_filename); 
end

mkdir results
save("results/" + source + "_" + "results2.mat", "eps_BP", "sigma_BP")
toc


% % Compare Sim to PA.

% pa = load("results/simulation_results.mat");
% sim = load("results/PA_results.mat");
% bias = (sim.eps_BP == 0) * 50;
% im_pa = imagesc(pa.eps_BP + bias);
% colorbar
% figure
% im_sim = imagesc(sim.eps_BP + bias);
% colorbar

