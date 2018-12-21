% by Zheng. Reconstruct depth from gradients

function [filledDepth] = reconstruct(depth, fillRegion, Dx, Dy, Lap)
    [Dxx, Dxy] = gradient(Dx);
    [Dyx, Dyy] = gradient(Dy);
    N = length(depth(fillRegion)); % Points to fill
    H = size(depth,1);
    W = size(depth,2);
    b = zeros(N,1);
    
    A_i = [];
    A_j = [];
    A_val = [];
    
    % Building fillregion mapping: 2d -> 1d
    index1D = zeros(H,W);
    index1D(fillRegion) = 1:N;
    % Building triplet list for sparse A 
    for i = 1:H
        for j = 1:W
            if fillRegion(i,j) == 0
                continue;
            end
            C_ij = 0;
            index = index1D(i,j);
            % Vertical
            if i == 1
                % D^(2,j) - D^(1,j) = Dy(1,j)
                C_ij = C_ij - 1;
                b(index) = Dy(i,j);
                if fillRegion(i+1,j) == 0
                    b(index) = b(index) - depth(i+1,j);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i+1,j);
                    A_val(end+1) = 1;
                end
            elseif i == H
                % - D^(i,j) + D^(i-1,j) = - Dy(H,j)
                C_ij = C_ij - 1;
                b(index) = - Dy(i,j);
                if fillRegion(i-1,j) == 0
                    b(index) = b(index) - depth(i-1,j);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i-1,j);
                    A_val(end+1) = 1;
                end
            else
                % D^(i-1,j) + D^(i+1,j) - 2D^(i,j) = Dyy(i,j)
                C_ij = C_ij - 2;
                b(index) = b(index) + Dyy(i,j);
                if fillRegion(i+1,j) == 0
                    b(index) = b(index) - depth(i+1,j);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i+1,j);
                    A_val(end+1) = 1;
                end
                if fillRegion(i-1,j) == 0
                    b(index) = b(index) - depth(i-1,j);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i-1,j);
                    A_val(end+1) = 1;
                end
            end
            % Horizontal
            if j == 1
                % D^(i,j+1) - D^(i,j) = Dx(i,j)
                C_ij = C_ij - 1;
                b(index) = Dx(i,j);
                if fillRegion(i,j+1) == 0
                    b(index) = b(index) - depth(i,j+1);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i,j+1);
                    A_val(end+1) = 1;
                end
            elseif j == W
                % - D^(i,j) + D^(i,j-1) = - Dx(i,j)
                C_ij = C_ij - 1;
                b(index) = - Dx(i,j);
                if fillRegion(i,j-1) == 0
                    b(index) = b(index) - depth(i,j-1);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i,j-1);
                    A_val(end+1) = 1;
                end
            else
                % D^(i,j-1) + D^(i,j+1) - 2D^(i,j) = Dxx(i,j)
                C_ij = C_ij - 2;
                b(index) = b(index) + Dxx(i,j);
                if fillRegion(i,j+1) == 0
                    b(index) = b(index) - depth(i,j+1);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i,j+1);
                    A_val(end+1) = 1;
                end
                if fillRegion(i,j-1) == 0
                    b(index) = b(index) - depth(i,j-1);
                else
                    % Add to triplet
                    A_i(end+1) = index;
                    A_j(end+1) = index1D(i,j-1);
                    A_val(end+1) = 1;
                end
            end
            A_i(end+1) = index;
            A_j(end+1) = index1D(i,j);
            A_val(end+1) = C_ij;
        end
    end
    % Build Sparse A from triplets
    A = sparse(A_i, A_j, A_val, N, N);
    
    % Solve the sparse matrix
    estimatedDepth = A\b;
    filledDepth = depth;
    filledDepth(fillRegion) = estimatedDepth;
end