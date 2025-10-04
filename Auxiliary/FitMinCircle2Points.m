function [Rmin, Cmin, X] = FitMinCircle2Points(X)
% Fit a minimum enclosing circle to a set of 2 or at most 3 points in
% 2-space.
%
% INPUT
%   - X     : M-by-2 array of point coordinates, where M<=3.
%
% OUTPUT
%   - Rmin  : radius of the minimum enclosing circle. 
%   - C     : coordinates of the circle centroid. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%

N = size(X,1);
if N>3
    error("nput must a N-by-2 array of point coordinates, with N<=3")
end

%% Empty set
if isempty(X)
    Cmin = nan(1,2);
    Rmin = nan; 
    return
end

%% A single point
if N==1
    Cmin = X;
    Rmin = 0;
    return
end

%% Line segment
if N==2
    Cmin = (X(2,:) + X(1,:))/2;
    Rmin = norm(X(2,:) - X(1,:))/2;
    return
end

%% Remove duplicate vertices, if there are any
interpointDistanceTol = 1E-12; 
DMat = zeros(N,N);
DMat(1,2) = norm(X(1,:) - X(2,:)); DMat(2,1) = DMat(1,2);
DMat(1,3) = norm(X(1,:) - X(3,:)); DMat(3,1) = DMat(1,3);
DMat(2,3) = norm(X(2,:) - X(3,:)); DMat(3,2) = DMat(2,3);
isDuplicatePair = DMat<=interpointDistanceTol;
isDuplicatePair(1:(N+1):end) = false;
if any(isDuplicatePair(:))
    for i = 1:N
        
        if size(X,1)<=i, break; end

        idx = isDuplicatePair(i,:);
        idx(1:i) = false;

        isDuplicatePair(idx,:) = [];
        isDuplicatePair(:,idx) = [];

        DMat(idx,:) = [];
        DMat(:,idx) = [];
        
        X(idx,:) = [];

    end
    [Rmin,Cmin,X] = FitCircle2Points(X);
    return
end
 
%% Longest edge
[~, maxPair] = max(DMat(:));
[a,b] = ind2sub([N N], maxPair);

%% Solve for minimum enclosing circle
c = 1:N;
c([a b]) = [];
Dca = (X(a,:) - X(c,:))/DMat(c,a);
Dcb = (X(b,:) - X(c,:))/DMat(c,b);
if (Dca*Dcb')<=eps % from inscribed angle theorem: iff angle ACB is >= 90 deg (i.e., dot prod <=0), C is inside/on the circle with diameter AB; where AB is the longest edge
    Cmin = (X(a,:) + X(b,:))/2;
    Rmin = norm(X(a,:) - X(b,:))/2;
    X = X([a b],:);
else
    A = X(1,:) - X(3,:);
    B = X(2,:) - X(3,:);

    A2 = sum(A.^2);
    B2 = sum(B.^2); 

    d = 2*(A(1)*B(2) - A(2)*B(1));
    Cmin = [B(2)*A2 - A(2)*B2, A(1)*B2 - B(1)*A2]/d;
    Rmin = norm(Cmin);
    Cmin = Cmin + X(3,:);
end
