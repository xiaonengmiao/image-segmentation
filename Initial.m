% Get Graph Laplacian for preparing 
% X - Original Labeled Image
% Y - Unlabeled Image
% Initializations on input image
% Input: 
%       I1: Original Image 1
%       I1_seg: Segmented Image 1
%       I2: Original Image 2
% Output:
%       c1: Initialization of first phase
%       c2: Initialization of second phase
%       c3: Initialization of third phase
%
% Matlab command window
% Implemented by Haohan Li, hlibb@connect.ust.hk
% The Hong Kong University of Science and Technology
% Oct 2015
% 
% all rights reserved
function [c1,c2,c3] = Initial(I1,I1_seg,I2)
[m1, n1,~] = size(I1);
[m2, n2,~] = size(I2);
c1 = zeros(m1*n1+m2*n2,1);
c2 = zeros(m1*n1+m2*n2,1);
c3 = zeros(m1*n1+m2*n2,1);
for i = 1:m1
    for j = 1:n1
        if sum(I1_seg(i,j,:)) > 600
            c1((i-1)*n2+j,1) = 1;
        elseif sum(I1_seg(i,j,:)) > 100
            c2((i-1)*n2+j,1) = 1;
        else
            c3((i-1)*n2+j,1) = 1;
        end
    end
end
for i = 1:m1
    for j = 1:n1
        I1_seg(i,j,:) = [255*c1((i-1)*n2+j), 125*c2((i-1)*n2+j), c3((i-1)*n2+j)];
    end
end
figure, imshow(I1_seg), title('segmented image 1');