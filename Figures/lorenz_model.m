function [t, x, y, z] = lorenz_model( z0, sigma, rho, beta, tmax )
%LORENZ MODEL Returns time (t) and time series (x, y, z) of the Lorenz model.
%
%   The Lorenz model is a three dimensional dynamical system described by
%   variables x, y and z that are solutions to the coupled differential
%   equations:
%
%       dx/dt = sigma*(y-x)
%       dy/dt = x*(rho - z) - y
%       dz/dt = x*y - beta*z
%
%   By introducing the variables:
%
%       z1 = x
%       z2 = y
%       z3 = z
%
%   the equations can then be rewritten as:
%
%       dz1/dt = sigma*(z2-z1)
%       dz2/dt = z1*(rho - z3) - z2
%       dz3/dt = z1*z2 - beta*z3
%
%   The input variable z0 contains the initial values for [z1 z3 z3].
%   The parameters are given by sigma, rho and beta.
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
%
[t,z] = ode45(@f,tspan,z0,options);
%
% Create timeseries variables for the output
%
x = timeseries(z(:,1),t);
y = timeseries(z(:,2),t);
z = timeseries(z(:,3),t);

    %
    % Nested function to provide derivatives for the coupled
    % equations.
    %
    function dzdt = f(t,z)
        % Derivative function. 
        % sigma, rho and beta are provided by the outer function.
        % dz1/dt = sigma*(z2-z1)
        % dz2/dt = z1*(rho - z3) - z2
        % dz3/dt = z1*z2 - beta*z3
          dzdt = [ sigma * (z(2) - z(1))
              z(1) * (rho - z(3)) - z(2)
              z(1) * z(2) - beta * z(3) ];
    end
    
    %
    % Nested function to provide the Jacobian, i.e., the partial
    % deriviatives J_{ij} = df_i/dz_j), where the function f is the rhs
    % of the system of differentia equations: dz_i/dt = f_i(z).
    %
    function dfdz = J(t,z)
      % Jacobian function.
      % sigma, rho and beta are provided by the outer function.
      dfdz = [
          -sigma       sigma     0
          rho-z(3)     -1        -z(1)
          z(2)         z(1)      -beta
          ];
   end
end

