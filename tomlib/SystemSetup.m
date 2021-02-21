
classdef SystemSetup < handle
    properties
        % We are interested in the following Physical Space:
        doi = [-100, 100, -115, 115];  % xmin xmax ymin ymax. (mm units)
        source_num = 16;  % Number of sources
        src_Tx_N = 16;  % number of transmitters
        src_Rx_N = 16;  % number of recievers.
        % The center of the array is thought of as the origin in above DOI.
        ant_xy = SystemSetup.readAntennaLocations();
        ant_ind;
        res;  % resolution in m
        mask;  % mask of DOI
        maskp;  % inverse of mask
        x_axis;
        y_axis;
        cell_x;
        cell_y;

        % Parameters as defined by Lei Guo (Standard Units). Those are the
        % ones that will be used by solvers.
        probes_Tx; % [Y X] locations of antennas.
        x_dash;  % meshing of Y
        y_dash;  % meshing of X
        Nx;  % number of cells in meshing x_dash
        Ny;  % number of cells in meshing y_dash
        total_n;  % total number of cells in the grid.
        axis_x;  % x grid
        axis_y;  % y grid.
        probes_Rx_coord;  % Coordinate Location of Rx after meshing
        probes_Rx_index;  % Index Locations of Rx after meshing.
        
        freq;  % Frequency in Hz
        bg;  % Background DP.
        cal1;  % Cal1 DP.
        cal2;  % Cal2 DP.
    end
    
    methods(Static)
        function xy = readAntennaLocations()
           pm = path_manager();
           ant_loc = readtable(pm.gdrive + "/emvision/Algorithm/toml_data/system/antenna_locations.csv");
           xy = [ant_loc.x, ant_loc.y];
        end
    end
    
    methods
        function obj = SystemSetup(resolution, background, cal1, cal2)  % Constructor.
            % Takes in on parameter: resolution (res)
            obj.res = resolution;
            obj.bg = background;
            obj.freq = obj.bg.freq;
            if nargin == 2
                cal1 = nan;
                cal2 = nan;
            elseif nargin == 3
                cal2 = nan;
            end
            obj.cal1 = cal1;
            obj.cal2 = cal2;
            obj = obj.meshDOI();
            obj.ant_ind = obj.getIndexFromLocation(obj.ant_xy);
            obj.mask = obj.getMask();
            obj.maskp = obj.mask == 0;
            
            obj.probes_Tx = fliplr(obj.ant_xy) * 1e-3;
            obj.x_dash = obj.y_axis * 1e-3;
            obj.y_dash = obj.x_axis * 1e-3;
            obj.Nx = length(obj.x_dash);
            obj.Ny = length(obj.y_dash);
            obj.total_n = obj.Nx * obj.Ny;
            [obj.axis_x, obj.axis_y] = ndgrid(obj.x_dash, obj.y_dash);
            
            %  Finding indices and coorindates of Rx's according to the
            %  meshing that happened:
            tmp_x = ones(obj.total_n, 1);
            tmp_y = ones(obj.total_n, 1);
            obj.probes_Rx_coord = zeros(obj.source_num, 2);
            obj.probes_Rx_index = zeros(obj.source_num, 2);

            for kk = 1 : obj.source_num
                dist = sqrt((obj.axis_x - obj.probes_Tx(kk, 1) .* tmp_x) .^ 2 + (obj.axis_y - obj.probes_Tx(kk, 2) .* tmp_y) .^ 2);
                [~, min_pos] = min(dist);

                [index_x, index_y] = ind2sub([obj.Nx, obj.Ny], min_pos);

                obj.probes_Rx_coord(kk, 1) = obj.axis_x(min_pos);
                obj.probes_Rx_coord(kk, 2) = obj.axis_y(min_pos);

                obj.probes_Rx_index(kk, 1) = index_x;
                obj.probes_Rx_index(kk, 2) = index_y;
            end

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
