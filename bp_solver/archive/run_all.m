
% Run the solver against all datasets.
clc
clear

pm = path_manager();
path1 = fullfile(pm.home, 'my_data', 'healthy47m/sparams');
path2 = fullfile(pm.home, 'my_data', 'target47m/sparams');


files1 = dir(path1);
files2 = dir(path2);
names = [struct2cell(files1), struct2cell(files2)];
files = {'sdgf', 'dfg'}; % random intialization
counter = 1;
for i=1: length(names)
    fname = string(cell2mat(names(1, i)));
    if (fname~= ".") && (fname ~= "..")
        files{counter} = strcat(cell2mat(names(2, i)), '\', cell2mat(names(1, i)));
        counter = counter + 1;
    end
end


N = length(files);
number_of_freq_samples = 13;    %%% Get the results under the frequency samples from 700 MHz to 900 MHz
resolution = '2mm';    %%% Use the resolution of 2mm

switch resolution
    case '2mm'
eps_BP = zeros(N, number_of_freq_samples, 112, 94);
sigma_BP = zeros(N, number_of_freq_samples, 112, 94);
    case '4mm'
eps_BP = zeros(N, number_of_freq_samples, 56, 47);
sigma_BP = zeros(N, number_of_freq_samples, 56, 47);
end

parallelize = false;

if parallelize == true
parpool(20)
parfor file=1:N
    for kk = 1 : number_of_freq_samples
        [eps_BP(file, kk, :, :), sigma_BP(file, kk, :, :)] = func_BP_CST(kk, files{file}, resolution);
    end
        fprintf('The %ith datapoint is finished...\n', file);    
end
delete(gcp('nocreate'))
else
for file=1:N
    for kk = 1 : number_of_freq_samples
        [eps_BP(file, kk, :, :), sigma_BP(file, kk, :, :)] = func_BP_CST(kk, files{file}, resolution);
    end
        fprintf('The %ith datapoint is finished...\n', file);    
end
end

