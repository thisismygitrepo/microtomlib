
classdef DOI
    properties
        % We are interested in the following Physical Space:
        doi = [-100, 100, -115, 115];  % xmin xmax ymin ymax
        ant_xy = DOI.readAntennaLocations();
        ant_ind;
        res;
        mask;
    end
    
    methods(Static)
        function xy = readAntennaLocations()
           pm = path_manager();
           ant_loc = readtable(pm.gdrive + "/emvision/Algorithm/toml_data/system/antenna_locations.csv");
           xy = [ant_loc.x, ant_loc.y];
        end
    end
    
    methods
        function obj = DOI(res)
            obj.res = res;
            obj.ant_ind = obj.getAntennaIndicesFromLocations();
            obj.mask = obj.getMask();
        end
        
        function [x, y] = meshDOI(obj)

            % now let us mesh it with a given resolution
            
            xx = obj.doi(1):obj.res:obj.doi(2) - obj.res;
            yy = obj.doi(3):obj.res:obj.doi(4);

            [x, y] = meshgrid(xx, yy);
            % At this point, 4 things happend:
            % 1- The physical x axis has become the columns of the image.
            % 2- The physical y axis has now become the rows of the image.
            % 3- The orientation of physical y axis has been flipped.
            %    i.e. it now increases by going down, instead of up.
            % 4- The origin of the y axis is flipped. It now starts at the top.
            % Summary:
               % In image space: origin is top left.
               % In Physics space: origin is bottom left. 
            % Thus, the physical origin, is now located at indices (1, 1)
            % of both x and y arrays returnd by this function.
            
            % Note: This is okay. Do not try to fix it.
            % At imshow time, you can pass argument 'lower' to make it show
            % in the correct orientation that matches physical space.
        end

        function mask = getMask(obj)
            % Takes in Resolution in mm units and returns mask image.

            [x, y] = obj.meshDOI();
            a = 217; %  Simplistic assumptions about the area being an ellipse.
            b = 174;  % these are the major and minor diagonals.

            mask = (y .^ 2 / (0.5 * a) .^ 2 + x .^ 2 / (0.5 * b) .^ 2) <= 1;
        end
        
        
        function ant_ind = getAntennaIndicesFromLocations(obj)
            % takes in an array xy of size N x 2 has physical locations.
            % returns an array of the same size has indices in the meshed
            % image of these physical locations.
            
            xx = obj.doi(1):obj.res:obj.doi(2) - obj.res;
            yy = obj.doi(3):obj.res:obj.doi(4);
            xy = obj.ant_xy;
            ant_ind = zeros(size(obj.ant_xy));
            
            for i=1:length(xy)
                x = xy(i, 1);
                q = find(xx > x);
                ant_ind(i, 2) = q(1) - 1;
                
                y = xy(i, 2);
                q = find(yy > y);
                ant_ind(i, 1) = q(1) - 1;
            end
        end
        
    end
    
end


