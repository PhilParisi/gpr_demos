function [k] = CalcKernel(x,y, hp)
    % Squared Exponential Kernel Function
    %   kij = exp(-abs(x - y)^2)
    L = hp.L;
    sigma = hp.sig;
    % a large sigma is a tigther fit - better at/between training data but
    % blows up beyond the training data
    qty = ((x-y)^2) / (2*L^2);
    k = (sigma^2)*exp(-qty);

end

% Great Resource
% http://evelinag.com/Ariadne/covarianceFunctions.html