classdef GeoBottom < LayeredFluidHalfspace
% GeoBottom    Expands layers with properties that vary with depth into a series of isotropic layers
%   GeoBottom is used to characterize the bottom layers with geoacoustic parameters and to calculate
%   an effective reflection coefficient and bottom loss.
%
%   Syntax:
%     GeoBottom(cpw, rhow, h, cp1, cp2, ap1, ap2, cs1, cs2, as1, as2, rho1, rho2)
%
%   Robert Chen
%   10 Feb 2015

  methods

    function obj = GeoBottom(cpw, rhow, h, cp1, cp2, ap1, ap2, cs1, cs2, as1, as2, rho1, rho2)
    % Constructor

      if ap1 == 0
        ap1 = zeros(size(cp1));
      end
      if ap2 == 0
        ap2 = zeros(size(cp2));
      end
      if cs1 == 0
        cs1 = zeros(size(cp1));
      end
      if cs2 == 0
        cs2 = zeros(size(cp2));
      end
      if as1 == 0
        as1 = zeros(size(cp1));
      end
      if as2 == 0
        as2 = zeros(size(cp2));
      end

      % Get indices of layers where any of the properties aren't constant across the layer
      idx = horzcat(cp1 ~= cp2, ap1 ~= ap2, cs1 ~= cs2, as1 ~= as2, rho1 ~= rho2);

      % Initialize
      hh = [];
      cp = [];
      ap = [];
      cs = [];
      as = [];
      rho = [];
      N = NaN(1, 5);

      % Loop through all specified layers
      for k = 1:length(h)

        if any(idx(k, :)) && k ~= length(h)

          param = find(idx(k, :));

          if ismember(1, param)
            N(1) = floor(abs(cp2(k) - cp1(k)) + 1);
          end

          if ismember(2, param)
            N(2) = 11;
          end

          if ismember(3, param)
            N(3) = floor((abs(cs2(k) - cs1(k)) + 1)/2);
          end

          if ismember(4, param)
            N(4) = 11;
          end

          if ismember(5, param)
            N(5) = 11;
          end

          Nmax = max(N);

          hh = vertcat(hh, diff(linspace(0, h(k), Nmax+1))');
          % hh = vertcat(hh, linspace(h(k)/Nmax, h(k)/Nmax, Nmax)');
          cp = vertcat(cp, linspace(cp1(k), cp2(k), Nmax)');
          ap = vertcat(ap, linspace(ap1(k), ap2(k), Nmax)');
          cs = vertcat(cs, linspace(cs1(k), cs2(k), Nmax)');
          as = vertcat(as, linspace(as1(k), as2(k), Nmax)');
          rho = vertcat(rho, linspace(rho1(k), rho2(k), Nmax)');

        else

          hh = vertcat(hh, h(k));
          cp = vertcat(cp, cp1(k));
          ap = vertcat(ap, ap1(k));
          cs = vertcat(cs, cs1(k));
          as = vertcat(as, as1(k));
          rho = vertcat(rho, rho1(k));

        end

      end

      obj = obj@LayeredFluidHalfspace(cpw, rhow, hh, cp, ap, cs, as, rho);

    end

    function varargout = plotBottomLoss(obj, frequencies, angles)
    % Calculate bottom loss and then plot the results

      if nargin < 3
        angles = 0:0.1:90;
      end

      BL = obj.getBottomLoss(frequencies, angles);
      [angles, frequencies] = meshgrid(angles, frequencies);
      colorVal = 0:14;
      colorMap = 'jet';
      figure('color', 'w');
      pcolor(angles, frequencies, BL');
      shading interp;
      hold on;
      contour(angles, frequencies, BL', colorVal(2:end-1), 'color', 'k');
      hold off;
      xlabel('Grazing Angle (deg)');
      ylabel('Frequency (Hz)');
      title('Bottom Reflection Loss');
      % EZColormap(colorVal, colorMap);
      colormap(colorMap);
      colorbar;
      xlabel(findobj(get(gcf, 'children'), 'Tag', 'Colorbar'), '(dB)');
      axis square;
      grid on;
      set(gca, 'layer', 'top');
      if nargout == 1
        varargout = BL;
      end

    end

  end

end
