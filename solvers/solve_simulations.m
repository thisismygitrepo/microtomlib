function [] = solve_simulations(idx)
% this functions recieves index, and sovles the corresponding simulation
% problem. the it saves the results.

database_path = "../Parameters_for_solvers";
load(database_path + "/paths_names.mat", 'paths', 'names')

name = strtrim(names(idx, :));
path = strtrim(paths(idx, :));
path = "../" + path;
cal1_path = "../cals/base_model_with_CalPhantom1.s16p";
cal2_path = "../cals/base_model_with_CalPhantom2.s16p";

disp(pwd)

solvers = ["CSI", "BIM", "MRCSI", "DBIM"];
for i=1:4
    asolver = solvers(i);
    [eps, sig] = run_multiFreq_solvers(asolver, path, cal1_path, cal2_path, database_path);
    save("../results/" + name + "_" + asolver + ".mat", "eps", "sig");
end

end