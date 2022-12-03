function [R,C] = ApproxMinBoundSphereND(X)
% Find a near-optimal bounding sphere for a set of N points in 
% D-dimensional space using Ritter's algorithm. In 3-space, the resulting
% sphere is approximately 5% bigger than the optimal minimum radius sphere.
%
% INPUT
%   - X     : N-by-D array of point coordinates, N is the number of point
%             samples and D>1 is their dimensionality.
%
% OUTPUT
%   - R     : radius of the approximate bounding sphere
%   - C     : centroid of the approximate bounding sphere
%
% REFERENCES:
% [1] Ritter, J. (1990), 'An efficient bounding sphere', in Graphics Gems,
%     A. Glassner, Ed. Academic Press, pp.301–303
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Basic error checking
if nargin<1 || ~ismatrix(X) || size(X,2)<2
    error('Invalid entry for 1st input argument (X)')
end

[N,D] = size(X);
[U,Xo] = deal([]);
if N==1
    R = 0; 
    C = X;
    return
elseif N>D && D<50  % standardize point cloud orientation
    Xo = mean(X,1);
    dX = bsxfun(@minus,X,Xo);
    [U,~,~] = svd(dX'*dX);
    X = dX*U; % (U'*X')'
end
    
% Convex hull of the point set
F = convhull(X);
F = unique(F(:));
Xch = X(F,:);

% Find points with the most extreme coordinates
[~,idx_min] = min(Xch,[],1);
[~,idx_max] = max(Xch,[],1);

Xmin = Xch(idx_min(:),:);
Xmax = Xch(idx_max(:),:);

% Compute distances between bounding-box vertices
X1 = permute(Xmin,[1 3 2]);
X2 = permute(Xmax,[3 1 2]);
D2 = bsxfun(@minus,X2,X1);
D2 = sum(D2.^2,3);

% Select point pair with the largest distance
[D2_max,idx] = max(D2(:));
[i,j] = ind2sub(size(D2),idx);
Xmin = Xmin(i,:);
Xmax = Xmax(j,:);

% Initial radius (squared) and centroid
R2 = D2_max/4;
R = sqrt(R2);
C = (Xmin + Xmax)/2;

% Loop through the M points, adjusting position and radius of the sphere
Nch = size(Xch,1);
for i = 1:Nch
    
    di2 = sum((Xch(i,:) - C).^2);
    if di2>R2
        di = sqrt(di2);
        R = (R + di)/2;
        R2 = R^2;
        dR = di - R;
        C = (R*C + dR*Xch(i,:))/di;
    end
    
end

% Undo rotation is centroid is required
if nargout>1 && ~isempty(U)
    C = C*(U') + Xo; % (U*C')' + Xo 
end
