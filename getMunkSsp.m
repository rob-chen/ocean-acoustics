function ssp = getMunkSsp(depth, dsca)
% getMunkSsp    Idealized deep water sound speed profile
%   The Munk profile is an idealized ocean sound speed profile that illustrates many features
%   typical of deep-water propagation.
%
%   Syntax:
%     ssp = getMunkSsp()
%     ssp = getMunkSsp(depth)
%     ssp = getMunkSsp(depth, dsca)
%
%   Inputs:
%     depth = Nx1 vector of water column depths in meters. Default is (0:100:5000)' m
%     dsca = Depth of minimum sound speed (deep sound channel axis) in meters. Default is 1300 m
%
%   Output:
%     ssp = Nx2 matrix of water column depths in meters and sound speeds in meters/second
%
%   Robert Chen
%   09 May 2014
  
  if nargin < 1
    depth = (0:100:5000)';
  else
    depth = depth(:);
  end
  
  if nargin < 2
    if depth(end) >= 3000
      dsca = 1300;
    else
      error('You must specify the depth of minimum sound speed, dsca!');
    end
  end
  
  z = 2 * (depth - dsca) / dsca;
  speed = 1500 * (1 + 0.00737*(z - 1 + exp(-z)));
  ssp = [depth, speed];

end
