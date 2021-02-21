classdef DD < handle
    %DD Dielectric Distribution Class
    %   Accepts (eps_r, sig) images and freq.
    
    properties
        eps_r;
        sig;
        setup;
        fig;
    end
    
    methods
        function obj = DD(eps_r, sig, setup)
            %DD Construct an instance of this class
            %   Detailed explanation goes here
            obj.eps_r = eps_r;
            
            if nargin == 1  % only eps is sent
                sig = false;  % define sigma
            end
            if nargin == 2 || nargin == 1
                setup = 0; % 
            end
            
            obj.sig = sig;
            obj.setup = setup;

        end
        
        function obj = imshow(obj, res)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if nargin == 1  % only obj but no res.
                res = 1;
            end
            
            if obj.setup == 0
                freq = "unknown";
            else
                freq = obj.setup.freq;
            end
            
            obj.fig = figure('Name',"Dielctric Proprties @" + string(freq) + " Hz");
            obj.fig.Position = [100, 100, 1500, 500];
            
            eps = imresize(obj.eps_r, res);
            eps = imgaussfilt(eps, 2);

            subplot(1, 2, 1)
            
            if obj.setup == 0
                imagesc(eps);
                axis off;
            else
                imagesc(obj.setup.axis_y(:), obj.setup.axis_x(:), eps)
            end
            
            axis image;  % Rectangular Window.
            view([0 -90]);  % Flip the rows to make it down to up like Y axis.
            colormap(jet);
            caxis([25 65]);
            colorbar;
            title("Relative Permitivity")
            
            sigma = imresize(obj.sig, res);
            sigma = imgaussfilt(sigma, 2);

            subplot(1, 2, 2)
            
            if obj.setup == 0
                imagesc(sigma);
                axis off;
            else
                imagesc(obj.setup.axis_y(:), obj.setup.axis_x(:), sigma)
            end
            
            view([0 -90]);  % Flip the rows to make it down to up like Y axis.
            axis image;
            colormap(jet);
            colorbar;
            caxis([0 0.5]);
            title("Conductivity")
            
            pause(0.001);
        end
        
        function show_probes(obj)
            
        end
        
        function show_mask(obj)
            subplot(1, 2, 1)
            imagesc(obj.setup.mask .* eps);
            view([0 -90]);  % Flip the rows to make it down to up like Y axis.

            subplot(1, 2, 2)
            imagesc(obj.setup.mask .* eps);
            view([0 -90]);  % Flip the rows to make it down to up like Y axis.

        end
    end
end

