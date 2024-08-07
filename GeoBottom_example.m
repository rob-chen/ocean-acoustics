%% GeoBottom Example
%
%   Reference:
%     Choi, J.W., Dahl, P.H., (2004). Mid-to-High-Frequency Bottom Loss in the East China Sea
%       IEEE Journal of Oceanic Engineering, Vol. 29. No. 4
%
%   Robert Chen
%   10 Feb 2015

% Lower water column sound speed and density
cpw = 1519;
rhow = 1000;

% Sediment layer and substrate properties
h = [0.9; inf];
cp1 = [1557; 1635];
cp2 = [1625; 1635];
ap1 = [0.25; 0.25]; % attenuation units are db/m/khz
ap2 = [0.25; 0.25]; % attenuation units are db/m/khz
cs1 = [0; 0];
cs2 = [0; 0];
as1 = [0; 0];
as2 = [0; 0];
rho1 = [2000; 2000];
rho2 = [2000; 2000];

% Create an instance of GeoBottom
obj = GeoBottom(cpw, rhow, h, cp1, cp2, ap1, ap2, cs1, cs2, as1, as2, rho1, rho2);
obj.attenUnits = 'db/m/khz';

% Plot the bottom loss as a function of frequency and grazing angle
obj.plotBottomLoss(0:10:2000, 0:90);

% Get bottom loss and effective reflection coefficient
[BL, R] = obj.getBottomLoss(0:10:2000, 0:90);
