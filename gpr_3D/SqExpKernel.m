function [k] = SqExpKernel(x1,x2, hp)
    % Squared Exponential Kernel Function
    % Euclidean Distance where applicable (dim > 1
    %   kij = exp(-abs(|x1 - x2|)^2)

    % for 2.5D case, x and y will be be 1x2 vector

    % Hyperparams
    % a large sigma is a tigther fit - better at/between training data but
    % blows up beyond the training data

    % Calculate Scalar Kernel Value
    dist = sqrt( (x1(1)-x2(1))^2 + (x1(2)-x2(2))^2 ); %scalar
    k = (hp.sigma^2) * exp( -(dist^2)/(2*hp.L^2) ); %scalar

end

% Great Resource
% http://evelinag.com/Ariadne/covarianceFunctions.html