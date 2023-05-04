function Y = spharm_real(X,L)
%SPHARM_REAL Real basis generator of spherical harmonics.
%   Y = spharm_real(X,L) computes the real bases of spherical harmonics 
%   of degree L and order M = -L, ..., L, evaluated for each element
%   of X.  X must be 3D coordinates on a sphere.
%
%   See also LEGENDRE.

%   Ilwoo Lyu, ilwoolyu@unist.ac.kr
%   Release: Apr 10, 2021

[phi,theta] = cart2sph(X(:,1),X(:,2),X(:,3));
theta = pi / 2 - theta;

m = [1:L]';
n = length(theta);
Pm = legendre(L,cos(theta),'sch');

lconstant = (2*L+1)/(4*pi);
precoeff = repmat(sqrt(lconstant),1,n);

% m < 0
lY = precoeff.*Pm(2:end,:).*sin(m*phi');
% m = 0
cY = sqrt(lconstant) * Pm(1,:);
% m > 0
rY = precoeff.*Pm(2:end,:).*cos(m*phi');

Y = [lY(end:-1:1,:); cY; rY]';
