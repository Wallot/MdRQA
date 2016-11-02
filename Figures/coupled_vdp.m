function [t, x, y] = coupled_vdp( z0, mu, e1, e2, tmax )
%COUPLED_VDP Returns time (t) and time series (x, y) of two coupled van der
%Pol oscillators.
%
%   The two coupled van der Pol (vdp) oscillators are described by
%   variables x and y that are solutions to the coupled differential
%   equations:
%
%       d2x/dt2 = mu*(1-x^2)*dx/dt - x + e1*(x-y)
%       d2y/dt2 = mu*(1-y^2)*dy/dt - y + e2*(y-x)
%
%   The coupling is asymmetric when e1 ~= e2, but the damping, given by mu
%   is the same for both vdp oscillators. By introducing the variables:
%
%       z1 = x
%       z2 = y
%       z3 = dx/dt
%       z4 = dy/dt
%
%   The equations can then be rewritten as:
%
%       dz1/dt = z3
%       dz2/dt = z4
%       dz3/dt = mu*(1-z1^2)*z3 - z1 + e1*(z1-z2)
%       dz4/dt = mu*(1-z2^2)*z4 - z2 + e2*(z2-z1)
%
%   The input variable z0 contains the initial values for [z1 z3 z3 z4].
%   The parameters are given by mu, e1, and e2.
%   The maximum time is tmax.
%
%   Dan Mønster, July 2016.

tspan = [0; tmax];
%
% The Jacobian is calculated analytically in the function J(t,z) provided
% below.
%
options = odeset('Jacobian',@J);
%
% The ode15s solver seems to be the best choice.
%
[t,z] = ode15s(@f,tspan,z0,options);
%
% Create timeseries variables for the output
%
x = timeseries(z(:,1),t);
y = timeseries(z(:,2),t);

    %
    % Nested function to provide derivatives for the coupled
    % equations.
    %
    function dzdt = f(t,z)
        % Derivative function. 
        % mu, e1, e2 is provided by the outer function.
        % dz1/dt = z3
        % dz2/dt = z4
        % dz3/dt = mu*(1-z1^2)*z3 - z1 + e1*(z1-z2)
        % dz4/dt = mu*(1-z2^2)*z4 - z2 + e2*(z2-z1)
          dzdt = [ z(3)
              z(4)
              mu*(1-z(1)^2)*z(3) - z(1) + e1*(z(1)-z(2))
              mu*(1-z(2)^2)*z(4) - z(2) + e2*(z(2)-z(1)) ];
    end
    
    %
    % Nested function to provide the Jacobian, i.e., the partial
    % deriviatives J_{ij} = df_i/dz_j), where the function f is the rhs
    % of the system of differentia equations: dz_i/dt = f_i(z).
    %
    function dfdz = J(t,z)
      % Jacobian function.
      % mu, e1, e2 is provided by the outer function.
      dfdz = [
          0                     0                       1               0
          0                     0                       0               1
          -2*mu*z(1)*z(3)-1+e1  -e1                     mu*(1-z(1)^2)   0
          -e2                   -2*mu*z(2)*z(4)-1+e2    0               mu*(1-z(2)^2)
          ];
   end
end

