function [Rmin,Cmin,Xb] = ExactMinBoundSphere3D(X)
% Compute exact minimum bounding sphere (MBS) of a 3D point cloud (or a 
% triangular surface mesh) using Welzl's algorithm [1]. 
%
% INPUT
%   - X     : M-by-3 list of point coordinates OR a triangular surface 
%             mesh specified as a (i) 'TriRep' object, (ii) 'triangulation'
%             object, (iii) 1-by-2 cell such that TR={Tri,V}, where Tri is
%             an M-by-3 array of faces and V is an N-by-3 array of vertex
%             coordinates, (iv) structure with 'faces' and 'vertices'
%             fields, such as the one returned by the 'isosurface'
%             function.
%
% OUTPUT
%   - Rmin  : radius of the minimum bounding sphere of X.
%   - Cmin  : 1-by-3 vector specifying centroid coordinates of the 
%             minimum bounding sphere of X.
%   - Xb    : K-by-3 array of MBS support points from which Rmin and Cmin
%             were derived; 2<=K<=4. Xb is a subset of X. See function 
%             'FitMinSphere2Points' for more info.
%
% REREFERENCES
% [1] Welzl, E. (1991), 'Smallest enclosing disks (balls and ellipsoids)',
%     Lecture Notes in Computer Science, Vol. 555, pp. 359-370
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if isnumeric(X) && ismatrix(X)
    if ~all(isfinite(X(:)))
        error("Point coordinates contain NaN or Inf entries. Remove them and try again.")
    elseif size(X,2)~=3
       error("This function works only for 3D data") 
    end    
else % if not a point cloud, the X may be a mesh
    try
        [~,X] = GetMeshData(X);
        if size(X,2)==2
            [Rmin,Cmin,Xb] = ExactMinBoundCircle(X);
            return
        end
    catch
        error("Invalid format for 1st input argument ('X')")
    end    
end


%% Pre-process the input
if size(X,1)>3
    
    % Check points for co-planarity 
    Xo = mean(X,1);
    dX = bsxfun(@minus,X,Xo);
    [U,D] = svd((dX'*dX)/size(X,1),0);
    D = diag(D);
    if D(3) < 1E-14
        [Rmin,Cmin,Xb] = ExactMinBoundCircle(dX*U(:,1:2));
        Cmin = (U(:,1:2)*Cmin(:))' + Xo;
        Xb = bsxfun(@plus,(U(:,1:2)*Xb')',Xo);
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


%% Get the MBS
if size(X,1)<4
    [Rmin,Cmin,Xb] = FitMinSphere2Points(X); 
    return
end

idxPermute = randperm(size(X,1));
X = X(idxPermute(:),:); % randomly permute the point set
if size(X,1)<1E3
    try %#ok<*TRYNC> 
        
        % Center and radius of the sphere
        [Rmin,Cmin] = B_MinSphere(X,[]);
        
        % Coordinates of the points used to compute parameters of the 
        % minimum bounding sphere
        D = sum(bsxfun(@minus,X,Cmin).^2,2);
        [D,idx] = sort(abs(D - Rmin^2));
        Xb = X(idx(1:4),:);
        D = D(1:4);
        Xb = Xb(D/Rmin*100<1E-3,:);
        [~,idx] = sort(Xb(:,1));
        Xb = Xb(idx,:);
        return
    end
end
    
% If we got to this point, then either size(X,1)>=1E3 or recursion depth 
% limit was reached. So need to break-up point-set into smaller subsets and
% then recombine the results.
M = size(X,1);
targetSubsetSize = min(floor(M/4),300);
res = mod(M,targetSubsetSize);
numSubsets = ceil(M/targetSubsetSize);  
numPtsPerSubset = targetSubsetSize*ones(1,numSubsets);
if res>0, numPtsPerSubset(end) = res; end
 
if res<0.25*targetSubsetSize && res>0
    numPtsPerSubset(numSubsets-1) = numPtsPerSubset(numSubsets-1) + numPtsPerSubset(numSubsets);
    numPtsPerSubset(numSubsets) = [];
    numSubsets = numSubsets - 1;
end

X = mat2cell(X,numPtsPerSubset,3);
Xb = [];
for i = 1:numSubsets
    
    % Center and radius of the sphere
    [Rmin,Cmin,Xi] = B_MinSphere([Xb; X{i}], []);    
    
    % 40 points closest to the sphere
    if i<numSubsets
        D = abs(sum(bsxfun(@minus,Xi,Cmin).^2,2) - Rmin^2);
    else
        D = abs(sqrt(sum(bsxfun(@minus,Xi,Cmin).^2,2)) - Rmin);
    end
    [D,idx] = sort(D);
    Xb = Xi(idx(1:min(40,numel(D))),:);
    
end

%% Points on the bounding sphere
D = D(1:min(40,numel(D)));
Xb = Xb(D/Rmin*100<1E-3,:);
if size(Xb,1)>4, Xb = Xb(1:4,:); end
[~,idx] = sort(Xb(:,1));
Xb = Xb(idx,:);


    function [R,C,P] = B_MinSphere(P,B)
    %   - P : N-by-3 list of points
    %   - B : K-by-3 list of MBS support points

        if size(B,1)==4 || isempty(P)
            [R,C] = FitMinSphere2Points(B); % fit sphere to boundary points (B)
            return
        end
        
        % Remove the last (i.e., end) point, p, from the P list
        P_new = P(1:end-1,:);
        p = P(end,:);
            
        % Check if p is on or inside the bounding sphere. If not, it must be
        % part of the new boundary.
        [R,C,P_new] = B_MinSphere(P_new,B);    
        isUpdateMBS = ~isfinite(R) || R<=eps || norm(p - C)>(R + eps);
        
        if isUpdateMBS
            B = [p; B];
            [R,C] = B_MinSphere(P_new,B);
            P = [p; P_new];
        end
        
    end

end

