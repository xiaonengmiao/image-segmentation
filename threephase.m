%======================================================================
% Data segmentation using diffuse interface model on graphs
% image segmentation
% 
% Implemented by Haohan Li, hlibb@connect.ust.hk
% The Hong Kong University of Science and Technology
% Oct 2015
% 
% all rights reserved
%
%----------------------------------------------------------------------
% Usage of Variables
% input:
%       I           = any gray/double/RGB input image
%       M           = no. of iteration
%       dt          = size of timestep
%       eps         = epsilon
%       c           = coefficient of convex splitting
% output:
%       I           = segmentation
%
%----------------------------------------------------------------------
% Description: Have not thought yet, to be expected.
%
%----------------------------------------------------------------------
% Please see the HELP file for details
%======================================================================

%%
%-- Initializations on input image I1 and I2
    % Read image I1, I1_seg and I2 as a matrix
    I1 = imread('118_1884.jpg'); 
    I1_seg = imread('118_1884_seg.jpg');
    I2 = imread('118_1888.jpg');
    [m1, n1, ~] = size(I1);
    [m1_seg, n1_seg, ~] = size(I1_seg);
    [m2, n2, ~] = size(I2);

%-- End Initializations on input image I1, I1_seg and I2
 
%%
%-- Core function
%-- Nystrom extention
    %[phi,O] = Nystrom(I1,I2);
    
%-- End Nystrom extention

%%
%-- Core function
    % Pre-condition
    M = 500;
    dt = 0.01;
    eps = 0.1;
    c = 21;
    [n,m] = size(phi);
    D = zeros(m,1);
    [c1_0,c2_0,c3_0] = Initial(I1,I1_seg,I2); 
    a1 = c1_0'*phi;
    b1 = c1_0'.^2*phi;
    d1 = c1_0'.^3*phi;
    a2 = c2_0'*phi;
    b2 = c2_0'.^2*phi;
    d2 = c2_0'.^3*phi;
    a3 = c3_0'*phi;
    b3 = c3_0'.^2*phi;
    d3 = c3_0'.^3*phi;
    e1 = [zeros(n1*m1,1)',zeros(n2*m2,1)'].*(c1_0')*phi;
    e2 = [zeros(n1*m1,1)',zeros(n2*m2,1)'].*(c2_0')*phi;
    e3 = [zeros(n1*m1,1)',zeros(n2*m2,1)'].*(c3_0')*phi;
    for i = 1:m
        D(i) = 1+dt*(eps*(3/4*(1-O(i,i))^2)+c);
    end
    %-- Main loop
    a1_0 = zeros(size(a1));
    a2_0 = zeros(size(a1));
    a3_0 = zeros(size(a1));
    ite = 0;
    while ite < M
        ite = ite+1;
        for i = 1:m;
            a1_0(i) = 1/D(i)*((1-8*dt/eps+c*dt)*a1(i)+24*dt/eps*b1(i)-16*dt/eps*d1(i)+4*dt/eps*a2(i)-12*dt/eps*b2(i)+8*dt/eps*d2(i)+4*dt/eps*a3(i)-12*dt/eps*b3(i)+8*dt/eps*d3(i)-dt*e1(i));
            a2_0(i) = 1/D(i)*((1-8*dt/eps+c*dt)*a2(i)+24*dt/eps*b2(i)-16*dt/eps*d2(i)+4*dt/eps*a1(i)-12*dt/eps*b1(i)+8*dt/eps*d1(i)+4*dt/eps*a3(i)-12*dt/eps*b3(i)+8*dt/eps*d3(i)-dt*e2(i));
            a3_0(i) = 1/D(i)*((1-8*dt/eps+c*dt)*a3(i)+24*dt/eps*b3(i)-16*dt/eps*d3(i)+4*dt/eps*a1(i)-12*dt/eps*b1(i)+8*dt/eps*d1(i)+4*dt/eps*a2(i)-12*dt/eps*b2(i)+8*dt/eps*d2(i)-dt*e3(i));
        end
        a1 = a1_0;
        a2 = a2_0;
        a3 = a3_0;
        c1 = phi*a1';
        c2 = phi*a2';
        c3 = phi*a3';
        a1 = c1'*phi;
        b1 = c1'.^2*phi;
        d1 = c1'.^3*phi;
        a2 = c2'*phi;
        b2 = c2'.^2*phi;
        d2 = c2'.^3*phi;
        a3 = c3'*phi;
        b3 = c3'.^2*phi;
        d3 = c3'.^3*phi;
        e1 = [ones(n1*m1,1)',zeros(n2*m2,1)'].*(c1'-c1_0')*phi;
        e2 = [ones(n1*m1,1)',zeros(n2*m2,1)'].*(c2'-c2_0')*phi;
        e3 = [ones(n1*m1,1)',zeros(n2*m2,1)'].*(c3'-c3_0')*phi;
    end
    save C1.mat c1
    save C2.mat c2
    save C3.mat c3
    
%-- End of Main loop

%-- plot
    for i = 1:m1
        for j = 1:n1
            I1(i,j,:) = [255*c1((i-1)*n2+j), 125*c2((i-1)*n2+j), c3((i-1)*n2+j)];
        end
    end
    for i = 1:m2
        for j = 1:n2
            I2(i,j,:) = [255*c1(m1*n1+(i-1)*n2+j), 125*c2(m1*n1+(i-1)*n2+j), c3(m1*n1+(i-1)*n2+j)];
        end
    end
    figure(2)
    imshow(I1)
    figure(3)
    imshow(I2)

%-- End of Core function
 
 
 
 
 