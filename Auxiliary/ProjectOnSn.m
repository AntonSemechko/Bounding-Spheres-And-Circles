function X = ProjectOnSn(X)
% Project a set of points (or mesh vertices) onto Sn.
%
% INPUT:
%   - X : N-by-d array of point coordinates, where N is the number of
%         points. Alternatively, X may be a surface mesh.
%   
% OUTPUT:
%   - X : N-by-d array of point coordinates such that norm(X(i,:))=1.
%         Alternatively, X may be a mesh whose vertices have beee projected 
%         onto a unit sphere.
%
% AUTHOR: Anton Semechko
%


if ismatrix(X) && isnumeric(X) % point cloud
    
    X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
    
else % mesh
    
    try
        [Tri,X,fmt] = GetMeshData(X);
    catch
        error('Invalid entry for 1st input argument (X)')
    end
    
    X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
    
    if fmt==1
        X = triangulation(Tri,X);
    elseif fmt==2
        X = TriRep(Tri,X); %#ok<*DTRIREP>
    elseif fmt==3
        X = {Tri X};
    else
        X = struct('faces',Tri,'vertices',X);
    end
    
end

