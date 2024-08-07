%% EZRay Example
%
%   Robert Chen
%   19 Aug 2015

%% Inputs

% Munk idealized ocean sound speed profile (m/s)
ssp.range = [0; 5000; 10000];
ssp.value{1} = getMunkSsp(0:100:3000, 1300);
ssp.value{2} = getMunkSsp(0:100:3000, 900);
ssp.value{3} = getMunkSsp(0:100:3000, 500);
ssp.interpx = false;

% Bathymetry (m)
% Note: We're just creating a finer bathymetry vector for more detailed ray paths
bathymetry = [(0:100:10000)', interp1([0 10000], [500, 3000], (0:100:10000))'];

% Fixed point depth (m)
fixDepth = 100;

% Departure angles (deg)
angles = -89:89;

%% Compute acoustic ray paths

% Create an instance of EZPZRay (the interface to EZRay)
% Note: You can call EZRay directory, of course, but EZPZRay adds plotting functionality.
obj = EZPZRay(ssp, bathymetry, fixDepth, angles);

% Perform the ray trace
obj.trace();

% Plot noise notch angles
obj.plot(obj.notchAngles);

% Plot contour of sound speed profiles
obj.plotSsp();
