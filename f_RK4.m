function f = f_RK4(nSteps, h, f0, rhs)
% f_RK4 - Fourth-order Runge-Kutta integrator for system of ODEs
%   f = f_RK4(nSteps, h, f0, rhs)
%   nSteps: Number of steps
%   h     : Step size (scalar)
%   f0    : Initial state (row vector)
%   rhs   : Function handle defining the ODEs (rhs(f), where f is row vector)
%
% Returns an nSteps x N array, where N is the number of state variables.

    % Preallocate the solution array
    f = zeros(nSteps, numel(f0));
    f(1,:) = f0;

    for i = 2:nSteps
        k1 = rhs(f(i-1, :));
        k2 = rhs(f(i-1, :) + 0.5*h*k1);
        k3 = rhs(f(i-1, :) + 0.5*h*k2);
        k4 = rhs(f(i-1, :) + h*k3);

        f(i,:) = f(i-1,:) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end
