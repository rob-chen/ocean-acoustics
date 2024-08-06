classdef GeoSurface < handle
% GeoSurface    Application of the Kirchhoff surface wave scattering model
%
%   References:
%     Etter, Paul C., (2003). Underwater Acoustic Modeling and Simulation (3rd ed.).
%       New York, NY: Spon Press
%     Jones, Adrian D., et al., (2009). Modelling the acoustic reflection loss at the rough ocean
%       surface. Proceedings of Acoustics. Adelaide, Australia: Australian Acoustical Society
%
%   Robert Chen
%   10 Mar 2015

  properties

    cpw; % Water sound speed (m/s)
    windspeed; % Wind speed at 19.5 m above sea level (m/s)

  end

  properties (Dependent, SetAccess = protected)

    h_sig; % Significant wave height (m)
    h_rms; % Root-mean-square wave height (m)
    h_avg; % Average wave height (m)

  end

  methods

    function obj = GeoSurface(soundspeed, windspeed, windspeedUnits)
    % Constructor

      % Sound speed in water
      obj.cpw = soundspeed;

      % The default windspeed is 0 m/s
      if nargin < 2
        obj.windspeed = 0;
      else
        obj.windspeed = windspeed;
      end

      % If windspeed is given in knots, convert to m/s
      if nargin == 3
        if strcmpi(windspeedUnits, 'kt')
          obj.windspeed = obj.windspeed*1852/3600;
        end
      end

    end

    function [SL, mu] = getSurfaceLoss(obj, frequencies, angles)
    % Calculate the reflection loss and coefficient due to scattering from surface waves

      % Reshape input arguments
      A = deg2rad(angles(:));
      F = frequencies(:)';

      % Wavenumbers
      k = 2*pi*F/obj.cpw;

      % Reflection coefficients
      mu = obj.refcoef(k, A);

      % Reflection loss
      SL = -10*log10(abs(mu));

    end

    function mu = refcoef(obj, k, theta)
    % Kirchhoff model
    % Note: Assumes gaussian probability of surface elevations with variance, a^2
    % Relationship between rms surface roughness and rms wave height from Etter(2003)(pp.66)

      % Rms surface roughness
      a = obj.h_rms/2;

      % Rayleigh roughness parameter
      gamma = 2*a*bsxfun(@times, k, sin(theta));

      % Reflection coefficient
      mu = -exp(-gamma.^2/2);

    end

  end

  % Setters and getters...
  methods

    function h_sig = get.h_sig(obj)
    % Significant wave height (average of the one-third highest waves)
    % Etter(2003)(pp.35)

      % Windspeeds are assumed to be in m/s
      h_sig = 0.21*obj.windspeed^2/9.81;

    end

    function h_rms = get.h_rms(obj)
    % Root-mean-square wave height

      h_rms = 0.704*obj.h_sig;

    end

    function h_avg = get.h_avg(obj)
    % Average wave height

      h_avg = sqrt(pi)/2*obj.h_rms;

    end

  end

end
