classdef EZGeoBottom < GeoBottom
% EZGeoBottom    Extends the GeoBottom class by letting the user specify layers by material name
%
%   Syntax:
%     obj = EZGeoBottom(h, material)
%     obj = EZGeoBottom(h, material, cpw, rhow)
%
%   The input parameters are:
%     h = Nx1 array of layer thicknesses
%     material = Nx1 cell array of strings
%     cpw = compressional sound speed in the water layer at the interface (default is 1500 m/s)
%     rhow = mass density of water at the interface (default is 1000 kg/m^3)
%
%   Valid material names:
%     clay, silt, sand, gravel, moraine, chalk, limestone, basalt
%
%   Robert Chen
%   15 Feb 2015

  properties (SetAccess = protected)

    layerSummary;

  end

  methods

    function obj = EZGeoBottom(h, material, cpw, rhow)
      % Constructor

      % Make sure h and material are column arrays:
      h = h(:);
      material = (:);

      % Generalized material properties
      materials = {'clay', 'silt', 'sand', 'gravel', 'moraine', 'chalk', 'limestone', 'basalt'};
      clay = {1.00, 0.2, @(d)50, 1.0, 1.5};
      silt = {1.05, 1.0, @(d)80*d^0.3, 1.5, 1.7};
      sand = {1.10, 0.8, @(d)100*d^0.3, 2.5, 1.9};
      gravel = {1.20, 0.6, @(d)180*d^0.3, 1.5, 2.0};
      moraine = {1.30, 0.4, @(d)600, 1.0, 2.1};
      chalk = {1.60, 0.2, @(d)1000, 0.5, 2.2};
      limestone = {2.00, 0.1, @(d)1500, 0.2, 2.4};
      basalt = {3.50, 0.1, @(d)2500, 0.2, 2.7};

      % Set default values for density and compressional wave speed of water
      if nargin < 4, rhow = 1000; end
      if nargin < 3, cpw = 1500; end

      % Cumulative depths of the layers
      d = [0; cumsum(h)];

      % Initialize
      cp1 = NaN(length(h), 1);
      cp2 = NaN(length(h), 1);
      ap1 = NaN(length(h), 1);
      ap2 = NaN(length(h), 1);
      cs1 = NaN(length(h), 1);
      cs2 = NaN(length(h), 1);
      as1 = NaN(length(h), 1);
      as2 = NaN(length(h), 1);
      rho1 = NaN(length(h), 1);
      rho2 = NaN(length(h), 1);

      % Loop through all layers
      for k = 1:length(h)

        idx = strcmpi(materials, material{k});

        if any(idx)

          cp1(k) = eval([materials{idx} '{1}'])*1500;
          cp2(k) = cp1(k);
          ap1(k) = eval([materials{idx} '{2}']);
          ap2(k) = ap1(k);
          cs1(k) = eval([materials{idx} '{3}(d(k))']);
          cs2(k) = eval([materials{idx} '{3}(d(k+1))']);
          as1(k) = eval([materials{idx} '{4}']);
          as2(k) = as1(k);
          rho1(k) = eval([materials{idx} '{5}'])*1000;
          rho2(k) = rho1(k);

        else

          error([material{k} ' is not a recognized material!']);

        end

      end

      obj = obj@GeoBottom(cpw, rhow, h, cp1, cp2, ap1, ap2, cs1, cs2, as1, as2, rho1, rho2);

      obj.layerSummary = table(material, h, cp1, cp2, cs1, cs2, rho1, rho2, ap1, ap2, as1, as2);

    end

  end

end
