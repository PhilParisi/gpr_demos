function [K] = K_Function(A, B, hp, style)
% Run input vectors through the Kernel!


K = zeros(size(A,1),size(B,1)); %predefine K for code efficiency
    % size(X,1) gives rows, size(X,2) gives cols

% 2D
    % Inputs A is mx1 and B is nx1 (vectors)
    %   Output will be mxn

% 3D
    % Inputs A is mx2 and B is nx2 (vectors)
    %   Output will be mxn

for i = 1:length(A)             %Start with rows of A, i = 1, 2, ...
    for j = 1:length(B)         %Go thru cols of B, j = 1, 2, ...
        
        % Choose Kernel specified in argument
        if style == "exact"
            K(i,j) = SqExpKernel(A(i,:), B(j,:), hp);
        elseif style == "sparse"
            K(i,j) = SqExpKernelSparse(A(i,:), B(j,:), hp);
        else %exact
            K(i,j) = SqExpKernel(A(i,:), B(j,:), hp);
        end

    end
end


end