% Interpolates a demographic function and gives error stats
function [Nint, stats] = getInterp(tsky, Nsky, t, Nt)

% Assumptions and notes
% - Nt is the true demographic function at t
% - Nsky is some skyline at tsky points
% - interpolation is ZOH as piecewise-constant

% ZOH interpolate skyline to t grid
Nint = interp1(tsky, Nsky, t, 'previous');
trange = range(t);

% Errors on grid
e = Nt - Nint;
% Bias, mse and variance
stats.bias = trapz(t, e)/trange;
stats.mse = trapz(t, e.^2)/trange;
stats.var = stats.mse - stats.bias^2;