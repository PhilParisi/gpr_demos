function [k] = SqExpKernel(x,y,hp)
% Squared Exponential Kernel Function
% Inputs: two datapoints (x,y) and hyperparameters (hp)

% hyperparameters stored in structure
    % a large sigma is a tigther fit - better at/between training data but
    %   blows up beyond the training data

%%%%% 2D Squared Exponential
if size(x(1,:),2) == 1                              % data contain scalars

    % distance between two points
    dist = (x-y);                                   % scalar
    % kernel value
    k = (hp.sigma^2)*exp( -(dist^2)/(2*hp.L^2) );   % scalar

end

%%%%% 2D Squared Exponential
if size(x(1,:),2) == 2                              % data contain 1x2 vectors

    % distance between two points
    dist = sqrt( (x(1)-y(1))^2 + (x(2)-y(2))^2 );   % scalar
    % kernel value
    k = (hp.sigma^2)*exp( -(dist^2)/(2*hp.L^2) );   % scalar
    
end

% Great Resource
% http://evelinag.com/Ariadne/covarianceFunctions.html