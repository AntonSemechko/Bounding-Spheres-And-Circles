function [X,R] = ProjectOnSn(X)
% Project a set of points (or mesh vertices) onto Sn.
%
% INPUT
%   - X : N-by-d array of point coordinates, where N is the number of
%         points. Alternatively, X may be a surface mesh.
%   
% OUTPUT
%   - X : N-by-d array of point coordinates such that norm(X(i,:)) = 1.
%         Alternatively, X may be a mesh whose vertices have been projected 
%         onto a unit sphere.
%   - R : N-by-1 array magnitudes of the position vectors in X before
%         projection.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if ismatrix(X) && isnumeric(X) % X is a point cloud

    R = sqrt(sum(X.^2,2));
    X = bsxfun(@rdivide,X,R);
    
else % X is a mesh?
    
    try
        [Tri,X,fmt] = GetMeshData(X);
    catch
        error("Invalid entry for 1st input argument ('X')")
    end
    
    R = sqrt(sum(X.^2,2));
    X = bsxfun(@rdivide,X,R);

    X = FormatOutputMesh(Tri,X,fmt);
    
end
