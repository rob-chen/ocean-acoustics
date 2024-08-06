function [idx1, idx2, negSlopeChk, Ydelta, idxAlwaysAboveZero] = getZeroCrossings(Y)
% getZeroCrossings
%   Get indices corresponding to zero crossings of the columns of an n-dimensional matrix
%
%   Syntax:
%     [idx1, idx2, negSlopeChk, Ydelta, idxAlwaysAboveZero] = getZeroCrossings(Y)
%
%   Robert Chen
%   25 May 2012

  % Change all zero entries to a very small positive number
  idxExactlyZero = (Y == 0);
  Y(idxExactlyZero) = eps(1);

  % Find the indices of Y where a sign change occurs indicating a zero-crossing
  idx1 = xor((Y < 0), circshift((Y < 0), -1));
  idx2 = xor((Y > 0), circshift((Y > 0),  1));

  % Correct for the circshift command at the end points
  colStr = repmat({':'}, [1, ndims(Y)-1]);
  idx1(end, colStr{:}) = 0;
  idx1(Y == 0) = 1;
  idx2(1, colStr{:}) = 0;
  idx2(Y == 0) = 1;

  % Get column indices for which Y is always greater than zero
  idxAlwaysAboveZero = all(Y > 0, 1);

  % Keep track of zero crossings where the slope of Y is negative. negSlopeChk is a vector the same
  % size as Ydelta that has 0s where the slope is positive and 1s where it is negative
  Ydelta = Y(idx1) - Y(idx2);
  negSlopeChk = (Ydelta >= 0);

end
