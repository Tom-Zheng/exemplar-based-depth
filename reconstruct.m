% by Zheng. Reconstruct depth from gradients

function [filledDepth] = reconstruct(depth, fillRegion, Dx, Dy)
    [Dxx, Dxy] = gradient(Dx);
    [Dyx, Dyy] = gradient(Dy);
    N = length(depth(fillRegion)); % Points to fill
    H = size(depth,1);
    W = size(depth,2);
    A = zeros(N);              % Building parameter matrix for Ax = b 
    b = zeros(N,1);
    
    % build A by point, then flattern it.
    A_temp = zeros(H,W,N);
    b_temp = zeros(H,W);
    
    for i = 1:H
        for j = 1:W
            if fillRegion(i,j) == 0
                continue;
            end
            b_temp(i,j) = Dxx(i,j) + Dyy(i,j); 
            row = zeros(H,W);
            row(i,j) = -4;
            if fillRegion(i+1,j) == 0
                b_temp(i,j) = b_temp(i,j) - depth(i+1,j);
            else
                row(i+1,j) = 1;
            end
            if fillRegion(i-1,j) == 0
                b_temp(i,j) = b_temp(i,j) - depth(i-1,j);
            else
                row(i-1,j) = 1;
            end
            if fillRegion(i,j+1) == 0
                b_temp(i,j) = b_temp(i,j) - depth(i,j+1);
            else
                row(i,j+1) = 1;
            end
            if fillRegion(i,j-1) == 0
                b_temp(i,j) = b_temp(i,j) - depth(i,j-1);
            else
                row(i,j-1) = 1;
            end
            A_temp(i,j,:) = row(fillRegion);
        end
    end
    % Mapping to rows of A.
    for i = 1:N
        column = A_temp(:,:,i);
        A(:,i) = column(fillRegion);
    end
    % Flatten to get b
    b = b_temp(fillRegion);
    estimatedDepth = A\b;
    filledDepth = depth;
    filledDepth(fillRegion) = estimatedDepth;
    
    % error = sum(sum(sqrt((depth - filledDepth).^2)));
end

% µ¥Ôª²âÊÔ1 Reconstruct Gaussian

% depth = gausswin(100)*(gausswin(100))';
% fillRegion = zeros(100);
% fillRegion(25:75,25:75) = 1;
% fillRegion = fillRegion == 1;
% indexmap = [];
% filled = reconstruct(depth, fillRegion, indexmap);
% error = sum(sum(sqrt((depth - filled)^2)));
% imagesc(filled);