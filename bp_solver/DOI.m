
classdef DOI
    properties
        % We are interested in the following Physical Space:
        doi = [-100, 100, -115, 115];  % xmin xmax ymin ymax. 
        source_num = 16;
        % The center of the array is thought of as the origin in above DOI.
        ant_xy = DOI.readAntennaLocations();
        ant_ind;
        res;
        mask;
        x_axis;
        y_axis;
        cell_x;
        cell_y;
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
            % Takes in on parameter: resolution (res)
            obj.res = res;
            obj = obj.meshDOI();
            obj.ant_ind = obj.getIndexFromLocation(obj.ant_xy);
            obj.mask = obj.getMask();
        end
        
        function obj = meshDOI(obj)

            % now let us mesh it with a given resolution
            
            obj.x_axis = obj.doi(1):obj.res:obj.doi(2) - obj.res;
            obj.y_axis = obj.doi(3):obj.res:obj.doi(4);

            [obj.cell_x, obj.cell_y] = meshgrid(obj.x_axis, obj.y_axis);
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
            a = 217; %  Simplistic assumptions about the area being an ellipse.
            b = 174;  % these are the major and minor diagonals.

            mask = (obj.cell_y .^ 2 / (0.5 * a) .^ 2 + obj.cell_x .^ 2 / (0.5 * b) .^ 2) <= 1;
        end
        
        function loc = getLocationFromIndex(obj, cell_idx)
           y = obj.doi(3) + (cell_idx(:, 1) - 1) * obj.res;  % get y from row
           x = obj.doi(1) + (cell_idx(:, 2) - 1) * obj.res;  % get x from col
           loc = [x, y];
        end
 
        function ind = getIndexFromLocation(obj, xy)
           ind = zeros(size(xy), 'uint32');
           tmp_x = obj.x_axis >= xy(:, 1);  % each row of xy has N column comparisons.
           tmp_y = obj.y_axis >= xy(:, 2);  % each row of xy has N column comparisons.
           % find function only works per vector.
           [num_points, ~] = size(xy);
           for i=1: num_points
               col = find(tmp_x(i, :), 1, 'first') - 1;  % find acts per column
               row = find(tmp_y(i, :), 1, 'first') - 1;
               ind(i, :) = [row, col];
           end
        end
        
    end
    
end

