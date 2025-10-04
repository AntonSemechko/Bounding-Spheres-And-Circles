function [Rmin,Cmin,Xb] = ExactMinBoundCircle(X)
% Compute exact minimum bounding circle (MBC) of a 2D point cloud using 
% Welzl's algorithm [1]. 
%
% INPUT
%   - X     : M-by-2 list of point coordinates, where M is the total
%             number of points.
%
% OUTPUT
%   - Rmin  : radius of the minimum bounding circle of X.
%   - Cmin  : 1-by-2 vector specifying centroid coordinates of the minimum
%             bounding circle of X.
%   - Xb    : K-by-2 array of MBC support points from which Rmin and Cmin
%             were derived; 2<=K<=3. Xb is a subset of X. See function
%             'FitMinCircle2Points' for more info.
%
% REREFERENCES:
% [1] Welzl, E. (1991), 'Smallest enclosing disks (balls and ellipsoids)',
%     Lecture Notes in Computer Science, Vol. 555, pp. 359-370
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<1 || isempty(X) || ~isnumeric(X) || ~ismatrix(X) 
    error("Unrecognized format for 1st input argument ('X')")
elseif size(X,2)~=2
    error("This function only works for 2D data")
elseif any(~isfinite(X(:)))
    error("Point data contains NaN and/or Inf entries. Remove them and try again.")    
end

%% Pre-process the input
if size(X,1)>2
    
    % Check points for collinearity
    Xo = mean(X,1);
    dX = bsxfun(@minus,X,Xo);
    [U,D] = svd((dX'*dX)/size(X,1),0);
    D = diag(D);
    if D(2)<1E-15
        dx = dX*U(:,1);
        [dx_min,id_min] = min(dx);
        [dx_max,id_max] = max(dx);
        Rmin = (dx_max - dx_min)/2;
        Cmin = U(:,1)'*(dx_max + dx_min)/2 + Xo;
        Xb = X([id_min; id_max],:);
        return
    end
    
    % Get convex hull of the point set
    F = convhull(X);
    F = unique(F(:));
    X = X(F,:);    
end

try
    X = uniquetol(X,eps,'ByRows',true); % remove duplicates
catch
    % older version of Matlab; 'uniquetol' is unavailable
end


%% Get the MBC
if size(X,1)<3
    [Rmin,Cmin,Xb] = FitMinCircle2Points(X); 
    return
end

idxPermute = randperm(size(X,1)); 
X = X(idxPermute(:),:); % randomly permute the point set
if size(X,1)<1E3
    try %#ok<TRYNC>
        % Center and radius of the circle
        [Rmin,Cmin] = B_MinCircle(X,[]);
        
        % MBC support points
        D = sum(bsxfun(@minus,X,Cmin).^2,2);
        [D,idx] = sort(abs(D-Rmin^2));
        Xb = X(idx(1:4),:);
        D = D(1:4);
        Xb = Xb(D<1E-6,:);
        [~,idx] = sort(Xb(:,1));
        Xb = Xb(idx,:);
        return
    end
end
    
% If we got to this point, then either size(X,1)>=1E3 or recursion depth 
% limit was reached. So need to break-up point-set into smaller subsets and 
% then recombine the results.
M = size(X,1);
targetSubsetSize = max(min(floor(M/4),300),3);
res = mod(M,targetSubsetSize);
numSubsets = ceil(M/targetSubsetSize);  
numPtsPerSubset = targetSubsetSize*ones(1,numSubsets);
if res>0, numPtsPerSubset(end) = res; end
 
if res<0.25*targetSubsetSize && res>0
    numPtsPerSubset(numSubsets-1) = numPtsPerSubset(numSubsets-1) + numPtsPerSubset(numSubsets);
    numPtsPerSubset(numSubsets) = [];
    numSubsets = numSubsets - 1;
end

X = mat2cell(X,numPtsPerSubset,2);
Xb = [];
for i = 1:numSubsets
    
    % Center and radius of the circle
    [Rmin,Cmin,Xi] = B_MinCircle([Xb;X{i}],[]);    
    
    % 40 points closest to the circle
    if i<numSubsets
        D = abs(sum(bsxfun(@minus,Xi,Cmin).^2,2) - Rmin^2);
    else
        D = abs(sqrt(sum(bsxfun(@minus,Xi,Cmin).^2,2)) - Rmin);
    end
    [D,idx] = sort(D);
    Xb = Xi(idx(1:min(40,numel(D))),:);
    
end

%% Points on the bounding circle
D = D(1:min(40,numel(D)));
Xb = Xb(D/Rmin*100<1E-3,:);
if size(Xb,1)>3, Xb = Xb(1:3,:); end
[~,idx] = sort(Xb(:,1));
Xb = Xb(idx,:);


    function [R,C,P] = B_MinCircle(P,B)
    %   - P : N-by-3 list of points
    %   - B : K-by-3 list of MBC support points
        
        if size(B,1)==3 || isempty(P)
            [R,C] = FitMinCircle2Points(B); % fit circle to boundary points
            return
        end
        
        % Remove the last (i.e., end) point, p, from the list
        P_new = P(1:end-1,:);
        p = P(end,:);
            
        % Check if p is on or inside the bounding circle. If not, it must be
        % part of the new boundary.
        [R,C,P_new] = B_MinCircle(P_new,B); 
        isUpdateMBC = ~isfinite(R) || R<=eps || norm(p - C)>(R + eps);
        
        if isUpdateMBC
            B = [p; B];
            [R,C] = B_MinCircle(P_new,B);
            P = [p; P_new];
        end
        
    end


end

