function [D, Gamma, MSD] = MSDcalc (xy, timeVec)

% Function to calculate diffusion coefficient, gamma and MSD of a 2D 
% diffusion track.
%
% Inputs:
%   -xy: table of xy coordinates, each row being 1 time point.  Must be in
%   desired units (um, nm, etc)
%   -timeVec: time (in seconds)
%
% Outputs:
%   -D = diffusion coefficient = slope(MSD)/4
%   -Gamma = curvature of the MSD graph, assuming MSD = 4Dt^(gamma)
%   -MSD = table of MSD values
%
%  Bryan Heit, May 2012

%% Check Inputs
if nargin ~= 2
    error ('Input is (xy, timeVec)');
end

if ndims(xy) ~= 2 || size(xy,2) ~= 2
    error ('Variable xy must be 2D array, 2 x n size');
end

MSD(1:(length(xy)-1)) = 0;

%% Generate MSD table
for i=1:(length(xy)-1)
    distTemp = (distyx(xy(i,:), xy(i+1:end,:))).^2; 
    MSD(1:length(distTemp)) = MSD(1:length(distTemp)) + distTemp;
    clear distTemp    
end

for i=1:(length(xy)-1)
    MSD(i) = MSD(i)/(length(xy)-i);
end

%% Calculate Gamma and D
f = ezfit ([timeVec:timeVec:(timeVec*(length(xy)-1))], MSD, 'y(x) = a*x^n');
Gamma = f.m(2);
D = f.m(1)/4;

end %end of function
