% Get Graph Laplacian for images
% Input: 
%       I1: Original Image 1
%       I2: Original Image 2
% Output:
%       phi: Eigenvectors of Graph Laplacian
%       O: Eigenvalues of Graph Laplacian
%
% Matlab command window
% Implemented by Haohan Li, hlibb@connect.ust.hk
% The Hong Kong University of Science and Technology
% Oct 2015
% 
% all rights reserved
function [phi,O] = Nystrom(I1,I2)
q = 0.1;
l = 20;
[m1, n1, ~] = size(I1);
[m2, n2, ~] = size(I2);
figure, imshow(I1), title('original image1');
figure, imshow(I2), title('original image2');
I1 = im2double(I1);
I2 = im2double(I2);
I1 = rgb2gray(I1);
I2 = rgb2gray(I2);
C = zeros(m1*n1+m2*n2,1);
for i = 1:m1
    for j = 1:n1
        C((i-1)*n1+j, 1) = I1(i, j);
    end
end
for i = 1:m2
    for j = 1:n2
        C((i-1)*n2+j+m1*n1, 1) = I2(i, j);
    end
end
X = random('unid',(m1*n1+m2*n2),l,1);
X = unique(X);
l = length(X);
Y = setdiff(1:(m1*n1+m2*n2),X);
Y = Y';
W11 = zeros(l,l);
W12 = zeros(l,(m1*n1+m2*n2-l));
for i = 1:l
    for j = 1:l
        if i == j
            W11(i,j) = 1;
        else
            W11(i,j) = exp(-norm([C(1+mod(X(i)-n1-1,m1*n1+m2*n2)), C(1+mod(X(i)-2,m1*n1+m2*n2)), C(1+mod(X(i)-1,m1*n1+m2*n2)), C(1+mod(X(i),m1*n1+m2*n2)), C(1+mod(X(i)+n1-1,m1*n1+m2*n2))]-[C(1+mod(X(j)-n1-1,m1*n1+m2*n2)), C(1+mod(X(j)-2,m1*n1+m2*n2)), C(1+mod(X(j)-1,m1*n1+m2*n2)), C(1+mod(X(j),m1*n1+m2*n2)), C(1+mod(X(j)+n1-1,m1*n1+m2*n2))])^2/q);
        end
    end
end
for i = 1:l
    for j = 1:(m1*n1+m2*n2-l)
        W12(i,j) = exp(-norm([C(1+mod(X(i)-n1-1,m1*n1+m2*n2)), C(1+mod(X(i)-2,m1*n1+m2*n2)), C(1+mod(X(i)-1,m1*n1+m2*n2)), C(1+mod(X(i),m1*n1+m2*n2)), C(1+mod(X(i)+n1-1,m1*n1+m2*n2))]-[C(1+mod(Y(j)-n1-1,m1*n1+m2*n2)), C(1+mod(Y(j)-2,m1*n1+m2*n2)), C(1+mod(Y(j)-1,m1*n1+m2*n2)), C(1+mod(Y(j),m1*n1+m2*n2)), C(1+mod(Y(j)+n1-1,m1*n1+m2*n2))])^2/q);
    end
end
dX = W11*ones(l,1)+W12*ones(2*m1*n1-l,1);
dY = W12'*ones(l,1)+W12'*(inv(W11))*(W12*ones(m1*n1+m2*n2-l,1));
W=(W11)^(-1);
sX = sqrt(dX);
sY = sqrt(dY);
W11 = W11./(sX*sX');
W12 = W12./(sX*sY');
[BX,T,CX] = svd(W11);
S = BX*(sqrt(inv(T)))*BX';
Q = W11+S*(W12*(W12)')*S;
[AX,O,DX] = svd(Q);
phi = [BX*sqrt(T);W12'*BX*sqrt(inv(T))]*(BX')*AX*(sqrt(inv(O)));
save phi.mat phi
save O.mat O