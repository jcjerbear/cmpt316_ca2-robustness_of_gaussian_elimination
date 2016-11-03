% MACM316 - Computing Assignment 2
% Title: Gaussian Elimination for a random matrix demo
% Description: Computes the mean solution error over Ntr trials for the
% system Ax=b where A is a random NxN matrix and x is a vector of ones.
% Plots the result as a histogram.

N = [10 13 15 17 21 27 29 33 39 41 47 52 55 62 69 73 80 91 96 103 112 129 136 140 150]; % Matrix size
Ntr = 300; % Number of trials

for i = 1:length(N)
    errs = zeros(Ntr,1); % 300x1 column vector of errors
    x = ones(N(i),1); % exact solution vector
    for j = 1:Ntr    
        A = randn(N(i),N(i)); % Construct a random NxN matrix (normally distributed)
        b = (A'*A)*x; % Compute the right-hand side vector
        z = (A'*A)\b; % Solve the linear system

        errs(j) = max(abs(z-x)); % Compute the error
    end
    % Compute the mean and standard deviation of the error
    mean_err(i) = mean(errs)
    sdev_err = sqrt(var(errs))
end

p = polyfit(log10(mean_err),log10(N),1)
r = polyval(p,log10(mean_err))

% Plot a graph of the errors
plot(log10(mean_err), log10(N), 'x');
hold on;
plot(log10(mean_err), r);
hold off;
title(['Mean Error VS N'])
xlabel('log_{10}(Mean Error)')
ylabel('log_{10}(N)')
grid on