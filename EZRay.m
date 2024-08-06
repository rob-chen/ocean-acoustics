function results = EZRay(ssp, bth, fixDepth, angles, rangeStart, rangeStop, maxBounces)
% EZRAY   Range-dependent ray-tracer
%   EZRAY calculates acoustic ray paths in environments where bathymetry and sound speed profiles
%   vary with horizontal range. Bathymetry and ssps are treated as piecewise-linear with undefined
%   gradients at their vertices. The water column is thus divided into layers wherein the sound
%   speed gradient is constant. This simplifies the problem quite a bit since ray paths then become
%   circle arcs within each layer.
%
%   Syntax:
%     results = EZRay(ssp, bathymetry, fixDepth)
%             = EZRay(..., angles)
%             = EZRay(..., angles, rangeStart)
%             = EZRay(..., angles, rangeStart, rangeStop)
%             = EZRay(..., angles, rangeStart, rangeStop, maxBounces)
%
%   Input:
%     Note: Be sure that you have defined a sound speed profile at range = 0.
%     ssp.range   => [Mx1]
%     ssp.value   => {Mx1} = {[Nx2]; [Nx2]; ...}
%     ssp.interpx => logical
%     bathymetry  => [Px2]
%     fixDepth    => scalar
%
%   Output:
%     results.sspx;
%     results.sspz;
%     results.sspv;
%     results.sspg;
%     results.bathymetry;
%     results.fixDepth;
%     results.angles;
%     results.rangeStart;
%     results.rangeStop;
%     results.maxBounces;
%     results.data;
%
%
%   Robert Chen
%   04 Aug 2015


    %% Check Inputs
        
    % Set default values:
    if nargin < 7 || isempty(maxBounces)
        maxBounces = 30;
    end
    if nargin < 6 || isempty(rangeStop)
        rangeStop = bth(end, 1);
    end
    if nargin < 5
        rangeStart = bth(1, 1);
    end
    if nargin < 4
        angles = -89.5:0.5:89.5;
    end

    % Make sure angles are converted to radians:
    angles = deg2rad(angles);

    % Set the default floating-point error tolerance:
    floatTol = 1e-6;

    % Make sure the fixed point depth doesn't exceed the bottom depth at the
    % starting range:
    Bchk = interp1q(bth(:, 1), bth(:, 2), rangeStart);
    if fixDepth > Bchk
        error(['The fixed point depth exceeds the bottom depth (' num2str(Bchk) ') at the starting range!']);
    end

    
    %% Process SSPs
    % Get all the unique depth values from all ssps and the max bottom depth.
    % These will be the common depth values that all ssps will be linearly
    % interpolated/extrapolated to:
    % sspz = ceil(max(bth(:, 2)));
    % for k = 1:length(ssp.value)
    %     sspz = unique([sspz; round(ssp.value{k}(:, 1))]);
    % end
    % EDIT: No longer round the depths when creating the vector of unique depths.
    sspz = max(bth(:, 2));
    for k = 1:length(ssp.value)
        sspz = unique([sspz; ssp.value{k}(:, 1)]);
    end
    
    % Linearly interpolate/extrapolate all ssps to the unique depths:
    sspv = NaN(length(sspz), length(ssp.value));
    for k = 1:length(ssp.value)
        sspv(:, k) = interp1(ssp.value{k}(:, 1), ssp.value{k}(:, 2), sspz, 'linear', 'extrap');
    end

    % Make sure there is a ssp at the final bathymetry range:
    sspx = ssp.range(:);
    if sspx(end) < bth(end, 1)
        sspx = vertcat(sspx, bth(end, 1));
        sspv = horzcat(sspv, sspv(:, end));
    end

    % If the ssp input struct has a field 'interpx' and its value is TRUE, then
    % the user wants the ssps interpolated at all ssp and bathymetry ranges.
    if isfield(ssp, 'interpx') && ssp.interpx == true
        
        % Interpolate all ssps to all ssp and bathymetry ranges:
        % sspv = interp1q(sspx, sspv', bth(:, 1))';
        % sspx = bth(:, 1);
        % EDIT: The vector of unique ranges is now combined from all ssps + 
        %       bathymetry ranges instead of just the bathymetry ranges.
        uniqRange = unique([sspx; bth(:, 1)]);
        sspv      = interp1q(sspx, sspv', uniqRange)';
        sspx      = uniqRange;
    
    % Otherwise, don't interpolate. Duplicate the ssp vectors across all ranges.
    else
        
        uniqRange = unique([sspx; bth(:, 1)]);
        sspvNew   = NaN(size(sspv, 1), length(uniqRange));
        [~, b]    = ismember(sspx, uniqRange);
        for k = 1:length(b)
            sspvNew(:, b(k):end) = repmat(sspv(:, k), [1, length(uniqRange)-b(k)+1]);
        end
        sspv = sspvNew;
        sspx = uniqRange;
        
    end

    % Calculate all gradients:
    sspg = bsxfun(@rdivide, diff(sspv), diff(sspz));
    sspg(abs(sspg) < floatTol) = 0;

    
    %% Trace
    % Pre-allocate output struct:
    results.sspx       = sspx;
    results.sspz       = sspz;
    results.sspv       = sspv;
    results.sspg       = sspg;
    results.bathymetry = bth;
    results.fixDepth   = fixDepth;
    results.angles     = angles(:);
    results.rangeStart = rangeStart;
    results.rangeStop  = rangeStop;
    results.maxBounces = maxBounces;
    results.data       = cell(length(angles), 1);

    % Store the ray termination code.
    % NOTE: 0 = Max range reached
    %       1 = Max bounce count reached
    %       2 = Ray started at the surface and was headed up
    %       3 = Ray stopped moving horizontally (probably because it was vertical)
    %       4 = Ray was vertical or reflected backwards off the bottom
    results.terminateCode = zeros(length(angles), 1);

    tic

    % Trace each ray one-by-one:
    for k = 1:length(angles)
        
        % Pre-allocate:
        x_pos      = NaN(1000, 1);          % Horizontal position
        z_pos      = NaN(1000, 1);          % Vertical position
        th_pos     = NaN(1000, 1);          % Ray angle
        rayLength  = NaN(1000, 1);          % Ray length
        bounceInfo = NaN(maxBounces, 6);    % Bounce info
        tau_pos    = NaN(1000, 1);          % Ray travel time

        % Define start values for depth and angle and range:
        x_pos(1)     = rangeStart;          % Current horizontal position
        z_pos(1)     = fixDepth;            % Current vertical position
        th_pos(1)    = angles(k);           % Current ray angle
        rayLength(1) = 0;                   % Current ray length
        tau_pos(1)   = 0;                   % Current travel time

        % Initialize various internal variables:
        R       = NaN;                      % Radius of curvature
        bounces = 0;                        % Number of bounces that have occurred
        x_idx   = 1;                        % Horizontal distance index
        S       = 0;                        % Starting ray length

        % Perform the trace while the range is less than the max range:
        while x_pos(x_idx) < rangeStop
            
            % Stop the trace if the max number of bounces (top+bottom) has been
            % exceeded:
            if bounces == maxBounces
                results.terminateCode(k) = 1;
                break
            end
            
            % If more space is needed for the ray position arrays, append
            % another pre-allocated block of space:
            if ~isnan(x_pos(end))
                x_pos     = vertcat(x_pos,     NaN(1000, 1));
                z_pos     = vertcat(z_pos,     NaN(1000, 1));
                th_pos    = vertcat(th_pos,    NaN(1000, 1));
                rayLength = vertcat(rayLength, NaN(1000, 1));
                tau_pos   = vertcat(tau_pos,   NaN(1000, 1));
            end
            
            % Initialize values for the current iteration using results from
            % the previous iteration:
            x0     = x_pos(x_idx);
            z0     = z_pos(x_idx);
            theta0 = th_pos(x_idx);
            
            % Get the indices that bracket which ssp is used at the current
            % horizontal position:
            idxS0 = find(x0 >= sspx, 1, 'last');
            idxS1 = idxS0 + 1;
            ssp   = [sspz, sspv(:, idxS0)];
            
            % Get the ssp indices that bracket the current depth:
            % NOTE: The sound speed profile is divided into layers of constant slope:
            idx0 = find(z0 >= ssp(:, 1), 1, 'last');
            idx1 = min(find(z0 <= ssp(:, 1)));
            
            % Get the bathymetry indices that bracket the current horizontal
            % position:
            idxB0    = find(x0 >= bth(:, 1), 1, 'last');
            idxB1    = idxB0 + 1;
            bthx_del = bth(idxB1, 1) - bth(idxB0, 1);
            bthz_del = bth(idxB1, 2) - bth(idxB0, 2);
            
            if idx0 == idx1             % The ray is at the boundary of two layers (i.e. a ssp data point)...
                
                % No need to interpolate; the speed is simply the current ssp
                % data point:
                c0 = ssp(idx0, 2);
                
                % Change the bracket indices according to the direction the ray
                % is headed.
                % NOTE: theta<0, Ray is headed UP
                %       theta>0, Ray is headed DOWN
                if theta0 < 0           % If the ray is heading up, it's entering the layer above the boundary...
                    if z0 == 0          %   unless it's already at the ocean surface...
                        results.terminateCode(k) = 2;
                        break           %   in which case, terminate the trace...
                    end
                    idx0 = idx0 - 1;    %   otherwise, make the first index the top of that layer.
                elseif theta0 > 0       % If the ray is heading down, it's entering the layer below the boundary...
                    idx1 = idx1 + 1;    %   make the second index the bottom of that layer.
                else                    % Otherwise, the ray is horizontal...
                    if R > 0            % If the radius of curvature was positive, the ray was curving down...
                        idx1 = idx1 + 1;%   so, make the second index the bottom of the layer below the boundary.
                    elseif R < 0        % If the radius of curvature was negative, the ray was curving up...
                        idx0 = idx0 - 1;%   so, make the first index the top of the layer above the boundary.
                    else                % Otherwise...
                        theta0 = 1e-6;  %   Assume the ray was headed down. (Need to revisit this!) Does this only occur when the ray starts out horizontal at a ssp vertex?
                        idx1 = idx1 + 1;
                    end
                end
                
            else                        % The ray is within a layer...
                
                % Linearly interpolate the sound speed:
                c0 = (z0 - ssp(idx0, 1)) / (ssp(idx1, 1) - ssp(idx0, 1)) * (ssp(idx1, 2) - ssp(idx0, 2)) + ssp(idx0, 2);
                
            end
            
            % Sound speed gradient at the current depth:
            % NOTE: g0<0, Downward refracting
            %       g0>0, Upward refracting
            g0 = sspg(idx0, idxS0);
            
            % If gradient is zero, then the current layer is an iso-speed
            % layer. Set the isoSpeed flag to the appropriate value:
            % NOTE: 0 => Local ssp gradient is NOT zero.
            %       1 => Local ssp gradient IS zero AND ray is NOT horizontal.
            %       2 => Local ssp gradient IS zero AND ray IS horizontal.
            if g0 ~= 0
                isoSpeed = 0;
            else
                isoSpeed = 1;
            end
            
            % Sound speed if the ray were to exit the current layer:
            if theta0 > 0               % If the ray is headed down...
                z1 = ssp(idx1, 1);      %   tentatively make the exit depth equal to the depth at the bottom of the layer.
                c1 = ssp(idx1, 2);      %   tentatively make the exit speed equal to the speed at the bottom of the layer.
            elseif theta0 < 0           % If the ray is headed up...
                z1 = ssp(idx0, 1);      %   tentatively make the exit depth equal to the depth at the top of the layer.
                c1 = ssp(idx0, 2);      %   tentatively make the exit speed equal to the speed at the top of the layer.
            else                        % Otherwise, the ray is horizontal, so...
                if g0 < 0               % If the ssp is downward refracting...
                    z1 = ssp(idx1, 1);  %   tentatively make the exit depth equal to the depth at the bottom of the layer.
                    c1 = ssp(idx1, 2);  %   tentatively make the exit speed equal to the speed at the bottom of the layer.
                elseif g0 > 0           % If the ssp is upward refracting...
                    z1 = ssp(idx0, 1);  %   tentatively make the exit depth equal to the depth at the top of the layer.
                    c1 = ssp(idx0, 2);  %   tentatively make the exit speed equal to the speed at the top of the layer.
                else                    % Otherwise, the ray is horizontal in an iso-speed layer...
                    z1 = z0;            %   so its depth will not change...
                    c1 = c0;            %   and neither will its speed change...
                    isoSpeed = 2;       %   and the isoSpeed flag should be set appropriately.
                end
            end
            
            % Calculate the tenative horizontal distance traveled:
            if isoSpeed == 0            % The local ssp is NOT constant...
                
                % Snell parameter:
                zeta = cos(theta0)/c0;
                
                % Radius of curvature:
                % NOTE: R<0, Ray is curving UP
                %       R>0, Ray is curving DOWN
                R = -1/zeta/g0;
                
                % Check if the ray turns completely within the current layer:
                if zeta^2*c1^2 < 1      % Ray does NOT turn...
                    
                    % Angle of the ray when it exits the layer:
                    % NOTE: The 'sign' function assumes that the ray will
                    % always refract in the same angle direction as the
                    % incoming angle direction. If the ray starts horizontal,
                    % set the exit angle sign based upon the sign of the radius
                    % of curvature:
                    if theta0 ~= 0
                        theta1 = sign(theta0)*acos(c1*cos(theta0)/c0);
                    else
                        theta1 = sign(R)*acos(c1*cos(theta0)/c0);
                    end
                    
                else                    % The ray DOES turn...
                    
                    % Angle of the ray when it turns:
                    theta1 = 0;
                    
                end
                
                % Tentative horizontal distance traveled while it was in the
                % layer:
                x_del = R*(sin(theta1) - sin(theta0));
                
            else
                
                % The ray angle remains unchanged:
                theta1 = theta0;
                
                if isoSpeed == 1        %   and the ray is NOT horizontal...
                    
                    % Calculate the horizontal distance traveled when it exits
                    % the layer:
                    x_del = abs((z1 - z0)/tan(theta0));
                    
                elseif isoSpeed == 2    %   or the ray IS horizontal...
                    
                    % Tentatively assume the ray remains horizontal for the
                    % remainder of the trace:
                    x_del = rangeStop - x0;
                    
                end
                
            end
            
            % Terminate the trace if the ray has stopped moving horizontally:
            if x_del == 0
                results.terminateCode(k) = 3;
                break
            end
            
            % Check if horizontal distance at which the ray exits the current
            % layer is greater than the distance to the next bathymetry point.
            % If so, set the ray horizontal position equal to this bathymetry
            % point horizontal position and calculate the new ray depth:
            if x0 + x_del > bth(idxB1, 1)
                x_del = bth(idxB1, 1) - x0;
            end
            
            % Check if the horizontal distance at which the ray exits the
            % current layer is greater than the distance to the next ssp range:
            if x0 + x_del > sspx(idxS1, 1)
                x_del = sspx(idxS1, 1) - x0;
            end
            
            % Calculate the new ray angle and vertical distance traveled:
            if isoSpeed == 0            % The local ssp is NOT constant...
                
                % Recalculate the exit angle:
                theta1 = asin(x_del/R + sin(theta0));
                
                % The vertical distance traveled while the ray was in the
                % layer:
                z_del = R*(cos(theta0) - cos(theta1));
                
                % Calculate the tentative path length:
                S = abs(R*(theta1 - theta0));
                
            elseif isoSpeed == 1        % The local ssp IS constant...
                
                % The vertical component of the straight line path:
                z_del = x_del*tan(theta0);
                
                % Calculate the tentative path length:
                S = sqrt(x_del^2 + z_del^2);
                
            elseif isoSpeed == 2        % The local ssp IS constant and the ray is horizontal...
                
                % The ray doesn't travel in the vertical direction at all:
                z_del = 0;
                
                % Calculate the tentative path length:
                S = x_del;
                
            end
            
            % Tentative new position of the ray for this iteration:
            x1 = x0 + x_del;
            z1 = z0 + z_del;
            
            % Now get the bottom depth at the ray's end horizontal position by
            % linearly interpolating the bathymetry values:
            bthz1 = (x1 - bth(idxB0, 1)) / (bth(idxB1, 1) - bth(idxB0, 1)) * (bth(idxB1, 2) - bth(idxB0, 2)) + bth(idxB0, 2);
            
            % Travel time:
            if isoSpeed == 0
                tau = abs(1/g0 * log(c1/c0 * (1+sin(theta0)) / (1+sin(theta1))));
            else
                tau = S/c0;
            end
            
            % Determine if the ray hits the surface or the bottom...
            % If it hits the surface (with a tolerance to account for floating
            % point errors)...
            if z1 < floatTol
                
                % Surface reflection only applies to non-horizontal rays:
                if isoSpeed < 2
                    
                    % Ocean surface is flat horizontal so the angle just
                    % changes sign:
                    theta1 = -theta1;
                    
                    % Update the bounce count:
                    bounces = bounces + 1;
                    
                    % Store surface bounce info:
                    % NOTE: [index surfBounce bottBounce angleAtContact angleAtReflect grazingAngle]
                    bounceInfo(bounces, :) = [x_idx+1 1 0 -theta1 theta1 abs(theta1)];
                    
                end
                
            % If the ray's depth at the end of this iteration is deeper than
            % the bathymetry point at that horizontal position, then the ray
            % has entered the bottom and its reflection must be accounted
            % for...
            elseif z1 > bthz1 - floatTol
                
                % Angle of bottom slope:
                alpha = atan(bthz_del/bthx_del);
                
                % Calculate where the ray hits the bottom using geometric
                % relationship between circular paths and straight line bottom
                % sections:
                if isoSpeed == 0
                    
                    % First, calculate location of the center of the circle
                    % that defines the rays path within the layer:
                    xc = x0 - R*sin(theta0);
                    zc = z0 + R*cos(theta0);
                    
                    % Then, determine the point of intersection between the ray
                    % path (circle) and the bottom (vector) by writing the
                    % latter as a pair of parametric equations with variable
                    % 't' and rearrange to use the quadratic equation:
                    A = bthx_del^2 + bthz_del^2;
                    B = 2*bthx_del*(bth(idxB0, 1) - xc) + 2*bthz_del*(bth(idxB0, 2) - zc);
                    C = (bth(idxB0, 1) - xc)^2 + (bth(idxB0, 2) - zc)^2 - R^2;
                    
                    % Solve for the parametric variable 't' and then solve for
                    % the intersection point (x1, z1):
                    t          = [-B, sqrt(B^2 - 4*A*C)]*[1, 1; 1, -1]*1/2/A;
                    xx         = bthx_del*t + bth(idxB0, 1);
                    zz         = bthz_del*t + bth(idxB0, 2);
                    [x1, idxt] = min(xx(xx>=x0));
                    z1         = zz(idxt);
                    
                    % Ray angle when it hits the bottom:
                    theta1b = asin((x1-x0)/R + sin(theta0));
                    
                    % Recalculate path length:
                    S = abs(R*(theta1 - theta0));
                    
                else                    % For the iso-speed case, it's the intersection of two lines...
                    
                    % Determine the point of intersection of the ray path
                    % (vector) and bottom section (vector) via parametric
                    % equations:
                    A  = [1 -1; tan(theta0) -tan(alpha)];
                    B  = [bth(idxB0, 1) - x0; bth(idxB0, 2) - z0];
                    t  = A\B;
                    x1 = bth(idxB0, 1) + t(2);
                    z1 = bth(idxB0, 2) + tan(alpha)*t(2);
                    
                    % Ray angle when it hits the bottom:
                    theta1b = theta0;
                    
                    % Recalculate path length:
                    S = sqrt(x_del^2 + z_del^2);
                    
                end
                
                % Now calculate the reflection angle...
                theta1 = 2*alpha - theta1b;
                
                % and the grazing angle:
                theta1g = abs(theta1 - alpha);
                
                % Update the bounce count:
                bounces = bounces + 1;
                
                % Store the bottom bounce info:
                % NOTE: [index surfBounce bottBounce angleAtContact angleAtReflect grazingAngle]
                bounceInfo(bounces, :) = [x_idx+1 0 1 theta1b theta1 theta1g];
                
                % Revised travel time:
                if isoSpeed == 0
                    tau = abs(1/g0 * log(c1/c0 * (1+sin(theta0))/(1+sin(theta1b))));
                else
                    tau = S/c0;
                end
                
            end
            
            % Adjust for propagating floating-point calculate errors:
            if ssp(idx1, 1) - z1 < floatTol
                z1 = ssp(idx1, 1);
            elseif z1 - ssp(idx0, 1) < floatTol
                z1 = ssp(idx0, 1);
            end
            
            % Horizontal position counter:
            x_idx = x_idx + 1;
            
            % Store the rays current position:
            x_pos(x_idx)     = x1;      % Horizontal position
            z_pos(x_idx)     = z1;      % Vertical position
            th_pos(x_idx)    = theta1;  % Ray angles (rad)
            rayLength(x_idx) = S;       % Ray length
            tau_pos(x_idx)   = tau;     % Ray travel time
            
            % Terminate the trace if the ray is vertical or reflects backwards:
            if abs(theta1) >= pi/2
                results.terminateCode(k) = 4;
                break
            end
            
        end
        
        % Store outputs:
        results.data{k}.x_pos      = x_pos(1:x_idx);
        results.data{k}.z_pos      = z_pos(1:x_idx);
        results.data{k}.theta      = th_pos(1:x_idx);
        results.data{k}.bounceInfo = bounceInfo(1:bounces, :);
        results.data{k}.rayLength  = cumsum(rayLength(1:x_idx));
        results.data{k}.tau_pos    = cumsum(tau_pos(1:x_idx));
        
    end

    toc

end
