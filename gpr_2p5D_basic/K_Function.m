function [K] = K_Function(A,B, hp)
% Run input vectors through the Kernel!

% This Function is designed for 2.5D GPR
% Inputs should be Nx2, ux2, Mx2 etc. (2 dimensional inputs)
% Output should be NxM = K_Function(Nx2,Mx2)

% Euclidean distance used

    K = zeros(length(A),length(B)); %predefine K for code efficiency

    for i = 1:length(A)             %Start with row i = 1
        for j = 1:length(B)         %Go thru cols j = 1, 2, 3
            K(i,j) = SqExpKernel(A(i,:),B(j,:), hp);
        end
    end

end