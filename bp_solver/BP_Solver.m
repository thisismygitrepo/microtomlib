
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call the func_BP_CST to get the eps_BP (reconstructed permittivity) and the sigma_BP (reconstructed conductivity)
%%% func_BP_CST has three parameters. The parameter 'freq_num' is the index of frequency sample
%%% Totally 13 frequency samples are supplied in func_BP_CST, they are 700MHz, 750MHz, 800MHz, 850MHz, 900MHz, 950MHz,
%%% 1000MHz, 1050MHz, 1100MHz, 1150MHz, 1200MHz, 1250MHz, and 1300MHz
%%% freq_num = 1 means the first frequency sample is selected, which is 700MHz
%%% freq_num = 2 means the second frequency sample is selected, which is 750MHz
%%% freq_num = 3 means the third frequency sample is selected, which is 800MHz, and so on...
%%% target_filename is the name of the S-parameter files. For example, target_filename = '101309_1.s16p' means we build the permittiivity
%%% and conductivity for the file "101309_1.s16p"
%%% To get the results from all the brain models, you need to write your own loop to visit all the S-parameter files 
%%% The parameter resolution specifies the resolution of the reconstructed images. It has to be a string variable
%%% Only two options are supplied for the parameter of resolution, they are '2mm' and '4mm'
%%% Therefore, use resolution = '2mm', or use resolution = '4mm'
%%% The folder "Parameters" includes all the variables that you need to use in func_BP_CST
%%% Please include the folder "Parameters" in your current working path before using this code
%%% I suggest to use the frequency samples of 700MHz ~ 950MHz. 
%%% I suggest do not go to the frequency samples higher than 950 MHz because based on the BP results, they look not very good

%%% Codes are written by Dr. Lei Guo from the University of Queensland
%%% Email: l.guo3@uq.edu.au
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

home = java.lang.System.getProperty('user.home');
path1 = fullfile(home, 'my_data', 'healthy47m/sparams');
path2 = fullfile(home, 'my_data', 'target47m/sparams');


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

