
clc
clear

% Select interpretor before Python is loaded up.
pyenv("Version", "C:\Users\s4551072\.conda\envs\gpuenv\python.exe")
% status is still unloaded so far. Once you run py.<command> it is loaded.

% matlab automatically loads up modules required to run any command.
% Therefore, one must put evrything in path in advance so that they're
% loaded up when needed.

insert(py.sys.path,int32(0), '')  % inserts current path @ location [0] in the list of paths (first entry) of python system.
% This is so that when we ask Python to look for something it knows we're
% doing this from here.

list = py.list({'Handerson', 'haha'});
py.demo_module.search(list)  % no need to import the module, this is done by matlab.
