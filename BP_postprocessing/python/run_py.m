
clc
clear

% Select interpretor before Python is loaded up.
pyenv("Version", 'C:\\Users\\s4551072\\.conda\\envs\\gpuenv\\python.exe')
% status is still unloaded so far. Once you run py.<command> it is loaded.

% matlab automatically loads up modules required to run any command.
% Therefore, one must put evrything in path in advance so that they're
% loaded up when needed.
insert(py.sys.path,int32(0), '')  % inserts current path @ loc[0]

py.search('nwq')