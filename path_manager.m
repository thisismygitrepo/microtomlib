classdef path_manager
    %PATH_MANAGER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        machine;
        gdrive;
        home;
    end
    
    methods
        function obj = path_manager()
%             machine = computer;
%             if machine == "GLNXA64"
%                 obj.machine = "linux";
%             elseif machine == "MACI64"
%                 obj.machine = "mac";
%             else
%                 obj.machine = "windows";
%                 obj.gdrive = "G:\";
%             end
            obj.home = string(java.lang.System.getProperty('user.home'));
        end

        function outputArg = gethis(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

