function [Rmin,Cmin,X] = FitMinSphere2Points(X)
% Fit a minimum enclosing sphere to a set of 2, 3, or at most 4 points
% in 3-space. 
%
% INPUT
%   - X     : N-by-3 array of point coordinates, where 2<=N<=4.
%
% OUTPUT
%   - Rmin  : radius of the minimum enclosing sphere
%   - Cmin  : 1-by-3 vector containing coordinates of sphere's centroid. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


N = size(X,1);
if N>4 
    error("Input must a N-by-3 array of point coordinates; 2<=N<=4")
end

%% Empty set
if isempty(X)
    Cmin = nan(1,3);
    Rmin = nan; 
    return
end

%% Single point
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
if N==4
    DMat(1,4) = norm(X(1,:) - X(4,:)); DMat(4,1) = DMat(1,4);
    DMat(2,4) = norm(X(2,:) - X(4,:)); DMat(4,2) = DMat(2,4);
    DMat(3,4) = norm(X(3,:) - X(4,:)); DMat(4,3) = DMat(3,4);
end
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
    [Rmin,Cmin,X] = FitSphere2Points(X);
    return
end

%% Three unique, though possibly collinear points
[~, furthestPointPair] = max(DMat(:));  % longest edge
[a,b] = ind2sub([N N], furthestPointPair);
cd = 1:N;
cd([a b]) = []; % remaining points

Dca = (X(a,:) - X(cd(1),:))/DMat(cd(1),a);
Dcb = (X(b,:) - X(cd(1),:))/DMat(cd(1),b);
if N==3
    if (Dca*Dcb')<=eps % from inscribed angle theorem: iff angle ACB is >= 90 deg (i.e., dot prod <=0), C is inside/on the circle with diameter AB; where AB is the longest edge
        Cmin = (X(a,:) + X(b,:))/2;
        Rmin = norm(X(a,:) - X(b,:))/2;
        X = X([a b],:);
    else
       [Rmin,Cmin] = getCircumCircle(X);
    end
    return
end


%% If we got to this point then we have 4 unique, though possibly collinear
% or coplanar points

% Check if the sphere that encloses the longest edge contains all points (rules out colinearity)
if (Dca*Dcb')<=eps
    Dda = (X(a,:) - X(cd(2),:))/DMat(cd(2),a);
    Ddb = (X(b,:) - X(cd(2),:))/DMat(cd(2),b);
    if (Dda*Ddb')<=eps
        Cmin = (X(a,:) + X(b,:))/2;
        Rmin = norm(X(a,:) - X(b,:))/2;
        X = X([a b],:);
        return
    end
end

% Check if the the points are coplanar
Dcd = (X(cd(2),:) - X(cd(1),:))/DMat(cd(1),cd(2));
vol = abs(cross(Dca,Dcb)*Dcd(:));
if vol < 1E-9
    
    % Look for circumcircle that contains all four points
    abcd = [a b cd([1 2]); ...
            a b cd([2 1]); ...
            cd([1 2]) a b; ...
            cd([1 2]) b a];

    [dR,Rmin_all] = deal(Inf(4,1));
    Cmin_all = nan(4,3);
    for i = 1:4
        abc = abcd(i,1:3);
        d = abcd(i,4);
        [Rmin_all(i),Cmin_all(i,:)] = getCircumCircle(X(abc,:));
        dR(i) = norm(Cmin_all(i,:) - X(d,:));
        if dR(i)<=(Rmin_all(i)+eps)
            X = X(abc,:);
            Rmin = Rmin_all(i);
            Cmin = Cmin_all(i,:);
            return
        end
    end

end


%% Ruled out collinearity and coplanarity, so minimum enclosing sphere must
% be supported by all four points 
AMat = 2*bsxfun(@minus,X(2:end,:),X(1,:));
b = sum(bsxfun(@minus,X(2:end,:).^2,X(1,:).^2),2);
Cmin = (AMat\b)';
Rmin = norm(X(1,:) - Cmin);


function [Rad,Cent] = getCircumCircle(X)
% Fit circumcircle to three points in 3-space 
%
%   - X     : 3-by-3 array of point coodinates, where X(i,:) corresponds to 
%             the i-th point

A = X(1,:) - X(3,:);
B = X(2,:) - X(3,:);
C = cross(A,B);

A2 = sum(A.^2);
B2 = sum(B.^2);
C2 = sum(C.^2);

Cent = cross(A2*B - B2*A,C)/(2*C2);
Rad  = norm(Cent);
Cent = Cent + X(3,:);

