classdef EZPZRay < handle
% EZPZRay   Interface to the range-dependent ray-tracer ezray.
%   EZPZRay wraps and extends the functionality of EZRay.
%
%   obj = EZPZRay(ssp, bathymetry, fixDepth)
%         EZPZRay(..., angles)
%         EZPZRay(..., angles, rangeStart)
%         EZPZRay(..., angles, rangeStart, rangeStop)
%         EZPZRay(..., angles, rangeStart, rangeStop, maxBounces)
%
%
%   Methods:
%
%      trace()
%
%      plot(angles)
%
%      plotNotch()
%
%      plotSsp()
%
%
%   Robert Chen
%   19 Aug 2015

    properties (SetObservable, AbortSet)
        
        angles;
        fixDepth;
        bathymetry;
        ssp;
        rangeStart;
        rangeStop;
        maxBounces;
        
    end
    
    properties (SetAccess = protected)
        
        pathInfo = {};
        bounceInfo = {};
        sspx;
        sspz;
        sspv;
        terminateCode;
        notchAngles;
        
    end
    
    methods

        function obj = EZPZRay(ssp, bathymetry, fixDepth, angles, rangeStart, rangeStop, maxBounces)
        % Constructor.
        
            % Store the input parameters:
            obj.ssp        = ssp;
            obj.bathymetry = bathymetry;
            obj.fixDepth   = fixDepth;
            
            % Set default ray angles:
            if nargin < 4
                obj.angles = -89.5:0.5:89.5;
            else
                obj.angles = angles;
            end
            
            % Set default start range:
            if nargin < 5
                obj.rangeStart = bathymetry(1, 1);
            else
                obj.rangeStart = rangeStart;
            end
            
            % Set default stop range:
            if nargin < 6 || isempty(rangeStop)
                obj.rangeStop = bathymetry(end, 1);
            else
                obj.rangeStop = rangeStop;
            end
            
            % Set default max number of bounces before the trace is stopped:
            if nargin < 7
                obj.maxBounces = [];
            else
                obj.maxBounces = maxBounces;
            end
            
            % Add listeners to reset the appropriate properties if any of the input parameters are changed:
            addlistener(obj, 'angles',     'PostSet', @obj.reset);
            addlistener(obj, 'fixDepth',   'PostSet', @obj.reset);
            addlistener(obj, 'bathymetry', 'PostSet', @obj.reset);
            addlistener(obj, 'ssp',        'PostSet', @obj.reset);
            addlistener(obj, 'rangeStart', 'PostSet', @obj.reset);
            addlistener(obj, 'rangeStop',  'PostSet', @obj.reset);
            addlistener(obj, 'maxBounces', 'PostSet', @obj.reset);
        
        end
        
        function trace(obj)
        % trace   Perform the ray trace using EZRay.

            % Call EZRay:
            disp('Tracing rays... please wait... ');
            rayInfo = EZRay(obj.ssp, obj.bathymetry, obj.fixDepth, obj.angles, obj.rangeStart, obj.rangeStop, obj.maxBounces);
            disp('Done!');

            % Store the output parameters:
            disp('Organizing the results... ');
            obj.sspx       = rayInfo.sspx;
            obj.sspz       = rayInfo.sspz;
            obj.sspv       = rayInfo.sspv;
            obj.pathInfo   = cell(length(obj.angles), 1);
            obj.bounceInfo = cell(length(obj.angles), 1);
            
            for k = 1:length(obj.angles)
                x                 = rayInfo.data{k}.x_pos;
                z                 = rayInfo.data{k}.z_pos;
                th                = rad2deg(rayInfo.data{k}.theta);
                rayLen            = rayInfo.data{k}.rayLength;
                travTime          = rayInfo.data{k}.tau_pos;
                obj.pathInfo{k}   = array2table(single([x, z, rayLen, travTime, th]), 'VariableNames', {'x', 'z', 'length', 'travelTime', 'angle'});
                obj.bounceInfo{k} = array2table(single([rayInfo.data{k}.bounceInfo(:, 1:3), rad2deg(rayInfo.data{k}.bounceInfo(:, 4:6))]), 'VariableNames', {'index', 'surfBounce', 'bottBounce', 'incidentAngle', 'reflectAngle', 'grazingAngle'});
            end
            
            % Store the termination codes:
            obj.terminateCode = rayInfo.terminateCode;
            
            % Get the rays that did not reach the surface:
            obj.notchAngles = obj.angles(cellfun(@(x) ~any(x.surfBounce), obj.bounceInfo) & obj.terminateCode < 2);
            
            disp('Done!');
            
        end
        
        function plot(obj, angles)
        % plot   Plot the ray paths.
        %
        %   obj.plot()
        %   obj.plot(angles);
            
            % By default, plot all rays:
            if nargin < 2
                angles = obj.angles;
            end
            
            % [idxRay, angles] = getNearest(obj.angles, angles, 0.0001);            
            % if all(isnan(idxRay))
            %     error('No matching angles found!');
            % else
            %     angles(isnan(angles)) = [];
            %     idxRay(isnan(idxRay)) = [];
            %     if length(angles) > 7
            %         Q = input(['You are about to display ' num2str(length(angles)) ' rays on a single plot. Are you sure? [N] '], 's');
            %         if ~strcmpi(Q, 'Y')
            %             return
            %         end
            %     end
            % end
            
            [idxRay, validAngles] = ismembertol(obj.angles, angles, 0.0001);
            angles = angles(validAngles(find(validAngles)));
            if ~any(idxRay)
                error('No matching angles found!');
            else
                idxRay = find(idxRay);
                % if length(angles) > 7
                %     Q = input(['You are about to display ' num2str(length(angles)) ' rays on a single plot. Are you sure? [N] '], 's');
                %     if ~strcmpi(Q, 'Y')
                %         return
                %     end
                % end
            end
            
            lineColors = lines(length(idxRay));
            
            figure('color', 'w', 'renderer', 'painters');
            
            h0 = patch([obj.bathymetry(:, 1); flipud(obj.bathymetry(:, 1))], ...
                       [obj.bathymetry(:, 2); repmat(1.1*max(obj.bathymetry(:, 2)), size(obj.bathymetry, 1), 1)], ...
                       [0.3, 0.3, 0.3]);
            set(h0, 'handlevisibility', 'off');
            set(gca, 'ydir', 'reverse');
            xlabel('Range (m)');
            ylabel('Depth (m)');
            hold on;
            
            h4 = NaN(length(idxRay), 1);
            
            for k = 1:length(idxRay)
                
                sbx = obj.pathInfo{idxRay(k)}.x(obj.bounceInfo{idxRay(k)}.index(obj.bounceInfo{idxRay(k)}.surfBounce == 1));
                sbz = obj.pathInfo{idxRay(k)}.z(obj.bounceInfo{idxRay(k)}.index(obj.bounceInfo{idxRay(k)}.surfBounce == 1));
                h1  = plot(sbx, sbz, 'ko', 'markerfacecolor', 'm', 'markersize', 4);
                set(h1, 'handlevisibility', 'off');
                
                bbx = obj.pathInfo{idxRay(k)}.x(obj.bounceInfo{idxRay(k)}.index(obj.bounceInfo{idxRay(k)}.bottBounce == 1));
                bbz = obj.pathInfo{idxRay(k)}.z(obj.bounceInfo{idxRay(k)}.index(obj.bounceInfo{idxRay(k)}.bottBounce == 1));
                h2  = plot(bbx, bbz, 'ko', 'markerfacecolor', 'y', 'markersize', 4);
                set(h2, 'handlevisibility', 'off');
                
                h4(k) = plot(obj.pathInfo{idxRay(k)}.x, obj.pathInfo{idxRay(k)}.z, 'o-', 'color', lineColors(k, :), 'markersize', 4);
                
            end
            
            for k = length(h4):-1:1
                uistack(h4(k), 'bottom');
            end
            
            legend(arrayfun(@(x) [num2str(x) '\circ'], angles, 'uni', 0), 'location', 'eastoutside');
            xlim([min(obj.bathymetry(:, 1)), max(obj.bathymetry(:, 1))]);
            ylim([0 1.1*max(obj.bathymetry(:, 2))]);
            set(gca, 'box', 'on', 'layer', 'top', 'color', [0.9, 0.99, 1], 'ticklength', [0 0]);
            grid on;
            hold off;
            
        end
        
        function plotNotch(obj)
        % plotNotch   Plot only non-surface-interacting rays.
        %
        %   obj.plotNotch()
            
            obj.plot(obj.notchAngles);
            
        end
        
        function H = plotSsp(obj, sspLims, sspLevels)
        % plotSsp   Plot the sound speed profiles.
        %
        %   obj.plotSsp(sspLims)
        %   obj.plotSsp(sspLims, sspLevels)
            
            figure('color', 'w', 'renderer', 'painters');
            
            % Define the default number of ssp colormap levels:
            if nargin < 3
                sspLevels = 31;
            end
            
            % Define the default ssp colormap limits:
            if nargin < 2
                sspLims = [floor(min(obj.sspv(:))), ceil(max(obj.sspv(:)))];
            end
            
            % Create the colormap:
            sspVals = linspace(sspLims(1), sspLims(2), sspLevels);
            colormap(jet(length(sspVals)-1));

            % Plot the sound speed profiles:
            [x, z] = meshgrid(obj.sspx, obj.sspz);
            h0     = pcolor(x, z, obj.sspv);
            set(h0, 'edgealpha', 0);
            
            % Colorbar:
            caxis([sspLims(1), sspLims(2)]);
            hcb = colorbar('ticks', sspVals);
            
            % Plot the bathymetry on top:
            hold on;
            h1 = patch([obj.bathymetry(:, 1); flipud(obj.bathymetry(:, 1))], ...
                       [obj.bathymetry(:, 2); repmat(1.1*max(obj.bathymetry(:, 2)), size(obj.bathymetry, 1), 1)], ...
                       [0.3, 0.3, 0.3]);
            set(h1, 'edgecolor', [0 0 0], 'zdata', ones(size(get(h1, 'xdata'))));
            ylim([0 1.1*max(obj.bathymetry(:, 2))]);
            hold off;
            
            % Finishing touches on the figure:
            xlabel('Range (m)');
            ylabel('Depth (m)');
            set(gca, 'ydir', 'reverse', 'layer', 'top', 'xgrid', 'on', 'ygrid', 'on', 'ticklength', [0 0]);
            
            % Return handles to graphics objects:
            H.h0  = h0;
            H.hcb = hcb;
            H.h1  = h1;
            
        end
    
    end
    
    methods (Access = protected)
        
        function reset(obj, src, evt)
        % Delete the old ray path and bounce data.
            
            obj.pathInfo      = {};
            obj.bounceInfo    = {};
            obj.sspx          = [];
            obj.sspz          = [];
            obj.sspv          = [];
            obj.terminateCode = [];
            obj.notchAngles   = [];
            
        end
        
    end
    
    methods
        
        function set.angles(obj, angles)
            
            obj.angles = angles(:);
            
        end
        
        
    end

end
