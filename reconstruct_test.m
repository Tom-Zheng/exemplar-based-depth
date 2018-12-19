depth = gausswin(100)*(gausswin(100))';
fillRegion = zeros(100);
fillRegion(50:100,50:100) = 1;
fillRegion = fillRegion == 1;

% [dx, dy] = gradient(depth);
% f = [0 1 0; 1 -4 1; 0 1 0];
% lap = imfilter(depth, f);

filled = reconstruct(depth, fillRegion, dx, dy);
error = sum(sum(sqrt((depth - filled)^2)));

imagesc(filled);

fprintf('Error = %f \n',error);