depth = gausswin(100)*(gausswin(100))';
fillRegion = zeros(100);
fillRegion(80:100,80:100) = 1;
fillRegion = fillRegion == 1;

[Dx,Dy] = gradient(depth);

filled = reconstruct(depth, fillRegion, Dx, Dy);

error = sum(sum(sqrt((depth - filled)^2)));
fprintf('Error : %f\n', error);
imagesc(filled);