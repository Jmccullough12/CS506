% CS-506 Homework3 - Problem 5
% Jeb and Jacob


original = imread('image.jpg');


[U1,S1,V1] = svdsketch(double(original),1e-1);
new1 = uint8(U1*S1*V1');


[U2,S2,V2] = svdsketch(double(original),1e-1,'MaxSubspaceDimension',45);
new2 = uint8(U2*S2*V2');


[U3,S3,V3] = svdsketch(double(original),1e-1,'MaxSubspaceDimension',35);
new3 = uint8(U3*S3*V3');

[U4,S4,V4] = svdsketch(double(original),1e-1,'MaxSubspaceDimension',25);
new4 = uint8(U4*S4*V4');

[U5,S5,V5] = svdsketch(double(original),1e-1,'MaxSubspaceDimension',15);
new5 = uint8(U5*S5*V5');



%new3 = uint8(U3*S3*V3');
%imshow(new3)
%title(sprintf('Rank %d approximation',size(S3,1)))

tiledlayout(2,3,'TileSpacing','Compact')
nexttile
imshow(original)
title(['Original (',sprintf('Rank %d)',rank(double(original)))])
nexttile
imshow(new1)
title(sprintf('Rank %d approximation',size(S1,1)))
nexttile
imshow(new2)
title(sprintf('Rank %d approximation',size(S2,1)))
nexttile
imshow(new3)
title(sprintf('Rank %d approximation',size(S3,1)))
nexttile
imshow(new4)
title(sprintf('Rank %d approximation',size(S4,1)))
nexttile
imshow(new5)
title(sprintf('Rank %d approximation',size(S5,1)))