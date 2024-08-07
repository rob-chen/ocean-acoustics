classdef LayeredFluidHalfspace < handle
% LayeredFluidHalfspace    Layered fluid halfspace model of ocean bottom layers
%   LayeredFluidHalfspace is used to characterize the bottom layers with geoacoustic parameters and
%   to calculate an effective reflection coefficient and bottom loss.
%
%   Syntax:
%     LayeredFluidHalfspace(cpw, rhow, h, cp, ap, cs, as, rho)
%
%   The following parameters are required to describe the water layer just above the ocean bottom:
%     cpw = compressional wave speed (m/s)
%     rhow = mass density (kg/m^3)
%
%   The following geoacoustic parameters are required of each bottom layer:
%     h = layer thickness (m)
%     cp = compressional wave speed (m/s)
%     ap = compressional wave attenuation (dB/lamda)
%     cs = shear wave speed (m/s)
%     as = shear wave attenuation (dB/lamda)
%     rho = mass density (kg/m^3)
%
%   Robert Chen
%   10 Feb 2015

  properties

    water; % Water properties
    bottom; % Bottom properties
    attenUnits; % Attenuation units

  end

  properties (Constant, Hidden)

    % Table variable names
    waterLayerVars = {'cp', 'rho'};
    bottomLayerVars = {'h', 'cp', 'ap', 'cs', 'as', 'rho'};

    % Valid attenuation coefficient units
    attenUnitsValid = {'db/lamda'; 'db/lambda'; 'db/m/khz'; 'db/m/hz'; 'db/km'; 'db/m'};

  end

  properties (Dependent, SetAccess = protected)

    numLayers; % Total number of bottom layers

  end

  methods

    function obj = LayeredFluidHalfspace(cpw, rhow, h, cp, ap, cs, as, rho)
    % Constructor

      % Store the specified water properties
      % obj.water = dataset([[cpw, rhow], obj.waterLayerVars]);
      obj.water = array2table([cpw, rhow], 'VariableNames', obj.waterLayerVars);

      % If any attenuation coefficients were defined to be zero, make them very small instead.
      % Otherwise, they can introduce instability to the reflection coefficient at low frequencies
      ap(ap == 0) = 1e-9;
      as(as == 0) = 1e-9;

      % Store the specified bottom properties
      % obj.bottom = dataset([[h, cp, ap, cs, as, rho], obj.bottomLayerVars]);
      obj.bottom = array2table([h, cp, ap, cs, as, rho], 'VariableNames', obj.bottomLayerVars);

      % Set the default attenuation units
      obj.attenUnits = 'db/lamda';

    end

    function [BL, Reff] = getBottomLoss(obj, frequencies, angles)
    % Calculate bottom loss and reflection coefficients

      % Reshape input arguments
      A = deg2rad(angles(:));
      F = frequencies(:)';

      % Compressional and shear wavenumbers
      kp = obj.wavenumber(F, reshape(obj.bottom.cp, 1, 1, obj.numLayers), ...
                             reshape(obj.bottom.ap, 1, 1, obj.numLayers), obj.attenUnits);
      ks = obj.wavenumber(F, reshape(obj.bottom.cs, 1, 1, obj.numLayers), ...
                             reshape(obj.bottom.as, 1, 1, obj.numLayers), obj.attenUnits);

      % Impedance of water layer above the ocean bottom
      z0 = obj.impedance(obj.water.rho, obj.water.cp, A);

      % The first angles to use are the specified angles of incidence at the water-bottom interface
      theta0p = A;
      theta0s = A;

      % The first wave numbers are those calculated from the compressional wave speed of water
      % (assuming no compressional attenuation)
      k0p = obj.wavenumber(F, obj.water.cp, 0, obj.attenUnits);
      k0s = k0p;

      % Pre-allocate the reflection coefficient and vertical phase delay matrices
      R = NaN(length(A), length(F), obj.numLayers-1);
      phi = NaN(length(A), length(F), obj.numLayers-1);

      % Progress through all the layers starting at the top, and calculate each layer's phase delay
      % and reflection coefficient with the layer above it
      for k = 1:obj.numLayers

        % Compressional wave refraction angles
        theta1p = obj.snell(theta0p, k0p, kp(:, :, k));

        % Compressional wave impedances
        z1p = obj.impedance(obj.bottom.rho(k), obj.bottom.cp(k), theta1p);

        % Effective impedance of the current layer
        % If the current layer doesn't support shear propagation, the effective impedance is just
        % the compressional impedance. Otherwise, the impedance is reduced, as shear propagation
        % acts as an additional degree of freedom for the acoustic wave to travel through the layer
        if obj.bottom.cs(k) == 0

          % When there is no shear, the effective impedance is the compressional impedance
          zeff = z1p;

        else

          % Shear wave refraction angles
          theta1s = obj.snell(theta0s, k0s, ks(:, :, k));

          % Shear wave impedances
          z1s = obj.impedance(obj.bottom.rho(k), obj.bottom.cs(k), theta1s);

          % Effective impedance
          zeff = z1p.*cos(2*theta1s).^2 + z1s.*sin(2*theta1s).^2;

        end

        % Reflection coefficient of the current layer
        R(:, :, k) = obj.refcoef(z0, zeff);

        % Check if the last layer has been reached
        if obj.numLayers == 1

          % If there is only one bottom layer (i.e., it's the substrate or halfspace)
          Reff = R;

        elseif k == obj.numLayers

          % If there is more than one layer and the current layer is the last one, work backwards to
          % calculate the effective reflection coefficient at the water-bottom interface
          Reff = R(:, :, k);
          for m = obj.numLayers:-1:2
            numer = R(:, :, m-1) + Reff.*exp(2i*phi(:, :, m-1));
            denom = 1 + R(:, :, m-1).*Reff.*exp(2i*phi(:, :, m-1));
            Reff = numer./denom;
          end

        else

          % Otherwise, calculate the vertical phase delay of the current layer
          phi(:, :, k) = obj.vertphasedelay(kp(:, :, k), obj.bottom.h(k), theta1p);

          % And prepare for the next iteration
          theta0p = theta1p;
          k0p = kp(:, :, k);
          z0 = zeff;

          if obj.bottom.cs(k) == 0
            theta0s = theta0p;
            k0s = k0p;
          else
            theta0s = theta1s;
            k0s = ks(:, :, k);
          end

        end

      end

      % Bottom loss (dB)
      BL = -20*log10(abs(Reff));

    end

  end

  methods (Static)

    function k = wavenumber(f, c, a, units)
    % Calculate wavenumbers (with attenuation): k = 2*pi*f/c + i*a
    % Note: Attenuation, a, must be in Np/m
    % The conversion factor is: 20*log10(e) = 8.686
    %   dB/lamda => a/c/8.686*f
    %   dB/m/kHz => a/1000/8.686*f
    %   dB/m/Hz => a/8.686*f
    %   dB/km => a/8.686/1000
    %   dB/m => a/8.686

      CF = 8.686;

      switch units

        case 'db/lamda'
          alpha = bsxfun(@times, a./c/CF, f);

        case 'db/lambda'
          alpha = bsxfun(@times, a./c/CF, f);

        case 'db/m/khz'
          alpha = bsxfun(@times, a/1000/CF, f);

        case 'db/m/hz'
          alpha = bsxfun(@times, a/CF, f);

        case 'db/km'
          alpha = a/1000/CF;

        case 'db/m'
          alpha = a/CF;

        otherwise
          error('Invalid attenuation units');

      end

      k = bsxfun(@plus, bsxfun(@rdivide, 2*pi*f, c), 1i*alpha);

    end

    function theta1 = snell(theta0, k0, k1)
    % Angle of refraction from Snell's Law: k0*cos(theta0) = k1*cos(theta1)

      theta1 = acos(bsxfun(@times, k0./k1, cos(theta0)));

    end

    function z = impedance(rho, c, theta)
    % Acoustic impedance: z = rho*c/sin(theta)

      z = rho.*c./sin(theta);

    end

    function R = refcoef(z0, z1)
    % Reflection coefficient between layer0 and layer1: R = (z1 - z0)/(z1 + z0)

      R = bsxfun(@minus, z1, z0)./bsxfun(@plus, z1, z0);

    end

    function phi = vertphasedelay(k, h, theta)
    % Vertical phase delay: phi = k*h*sin(theta)

      phi = h*bsxfun(@times, k, sin(theta));

    end

  end

  % Setters and getters...
  methods

    function numLayers = get.numLayers(obj)
    % Total number of defined bottom layers

      numLayers = size(obj.bottom, 1);

    end

    function set.attenUnits(obj, units)
    % Make sure the attenuation units are all valid and lowercase

      if ~any(strcmpi(obj.attenUnitsValid, units))

        disp('Invalid attenuation units!');
        disp('Valid units are:');
        disp(obj.attenUnitsValid);
        return

      end

      obj.attenUnits = lower(units);

    end

  end

end
