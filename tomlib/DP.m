classdef DP < handle  % see comparing handle class vs value class.
    % DP Dielctric (Electric) Properties class.
    %  Gives easy access to various forms of properties.
    
    properties
        path;
        freq;
        w;
        
        sig;  % conductivity
        eps;  % permitivity        
        eps_complex;  % complex permitivity
        
        background;  % background for relativity.
        eps_r;  % relative permitivity
        eps_complex_r;  % relative complex permitivity;
        contrast;  % contrast = eps_complex_r - 1.
        % To get any of below: use real(.) imag(.) with above
%         eps_p_r;  % real part of relative complex permitivity
%         eps_pp_r; % imaginary part of relative complex permitivity
%         eps_p; % real part of complex permitivity
%         eps_pp; % imaginary part of complex permitivity
            
        kb;  % propagation constant.
        eps_o = 8.854187817e-12;
        uo = 4e-7 * pi;  
        c = 299792458;
    end
    
% kb = w * sqrt(const.uo * eps_b);    end pm.join(pm.gdrive, "emvision\Algorithm\toml_data\system\Material Properties\
    
    methods
        function obj = DP(freq, path)
            %DP Construct an instance of this class
            %   Requires path to complex dispersive properties and frequency
            %   required.
            
            if nargin == 1 || nargin == 0
                % in case of no input, default to air as material
                obj.path = "air";
                
                if nargin == 0
                    obj.freq = 1;  % first input is frequency.
                else
                    obj.freq = freq;  % first input is frequency.
                end
                
                obj.w = 2 * pi * obj.freq;
                obj.sig = 0;
                obj.eps = obj.eps_o;
                obj.eps_complex = obj.eps - 1i * obj.sig / obj.w;
                
                obj.eps_r = 1;
                obj.eps_complex_r = obj.eps_r - 1j * obj.sig - 1;
                obj.contrast = obj.eps_complex_r - 1;
                obj.kb = obj.w * sqrt(obj.uo * obj.eps);
            else
            
            obj.path = path;
            obj.freq = freq;
            obj.w = 2 * pi * obj.freq;
            [eps_p_r, eps_pp_r] = getProperties(path, freq);
            % those are relative to air, thus, to get absolute values:
            
            obj.sig = eps_pp_r * obj.w * obj.eps_o;
            obj.eps = eps_p_r * obj.eps_o;
            
            % Computing relative values while assuming air as a background
            obj.relative_to(DP(obj.freq));
            
            end
            
        end
        
        function obj = relative_to(obj, background)
            %relative_to Populates relative properties of the material.
            %   requires another instance of the class.
            obj.eps_complex = obj.eps - 1i * obj.sig / obj.w;
            obj.kb = obj.w * sqrt(obj.uo * obj.eps_complex);
            
            obj.background = background;
            obj.eps_r = obj.eps / background.eps;
            obj.eps_complex_r = obj.eps_complex / background.eps_complex;
            obj.contrast = obj.eps_complex_r - 1;
           
        end
    end
    
    methods(Static)
        function obj = from_props(freq, eps, sig)
            obj = DP(freq);
            obj.path = "manual_insertion";
            obj.eps = eps;
            obj.sig = sig;
            obj.relative_to(DP(obj.freq));   
        end
    end
end

