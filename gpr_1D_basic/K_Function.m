function [K] = K_Function(A,B, hp)
% Run input vectors through the Kernel!
% Inputs should both be 1xm and 1xn arrays
%   Output will be mxn

    K = zeros(length(A),length(B)); %predefine K for code efficiency

    for i = 1:length(A)             %Start with row i = 1
        for j = 1:length(B)         %Go thru cols j = 1, 2, 3
            K(i,j) = CalcKernel(A(i),B(j), hp);
        end
    end

end