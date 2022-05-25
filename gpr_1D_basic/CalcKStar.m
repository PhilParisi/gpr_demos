function [K_Star] = CalcKStar(matX,newX)
% Calculate K_Star Matrix
%   K_Start output is vertical nx1

    for i = 1:length(matX)
        K_Star(i,1) = CalcKernel(matX(i),newX);
    end

end