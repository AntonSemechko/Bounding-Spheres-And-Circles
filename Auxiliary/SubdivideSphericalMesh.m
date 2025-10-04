function [TR, G, SM] = SubdivideSphericalMesh(TR, k, G, W)
% Subdivide triangular, quadrilateral, or mixed element surface mesh, 
% representing a zero-centered sphere, k times using triangular
% (or quadrilateral) quadrisection. See functions 'TriQuad' and 'QuadQuad'
% for more info.
%
% INPUT
%   - TR   : surface mesh of a ZERO-CENTERED sphere represented as an 
%            object of 'TriRep' class, 'triangulation' class, or a cell 
%            such that TR={F,X}, where F is an M-by-3 or M-by-4 array of
%            faces, and X is an N-by-3 array of vertex coordinates.
%   - k    : (optional) desired number of subdivisions. k=1 is default.
%   - G    : (optional) scalar, vector, or tensor field defined at the 
%            vertices of TR.
%   - W    : (optional) positive weights associated with vertices of TR.
%            See function 'TriQuad' (or 'QuadQuad') for more info.
%
% OUTPUT
%   - TR  : subdivided mesh. Same format as input mesh.
%   - G   : interpolated scalar, vector, or tensor field defined at the
%           vertices of the subdivided mesh. Note, linear interpolation of
%           tensor fields may produce output that is not a valid tensor.
%   - SM  : sparse K-by-N subdivision matrix, where K is the number of 
%           vertices in the subdivided mesh.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<2 || isempty(k), k = 1; end
if nargin<3, G = []; end
if nargin<4, W = []; end


% Get the data structure
[F,X,fmt] = GetMeshData(TR);
if fmt==1 && size(F,2)==4
    error('Tet meshes cannot be processed by this function')
end

% Make sure vertices of the input mesh lie on a unit sphere
R = sqrt(sum(X.^2,2));
X = bsxfun(@rdivide,X,R);

% Return mesh as is if k<1
k = round(k(1));
if k<1 
    TR = FormatOutputMesh(F,X,fmt);
    return
end

% Spherical subdivision
for i = 1:k
    
    % Subdivide the mesh
    if nargout<=2    
        if size(F,2)==3
            [TR,W,G] = TriQuad({F X},W,G);
        else
            [TR,W,G] = QuadQuad({F X},W,G);
        end
    else
        if size(F,2)==3
            [TR,W,G,SMi] = TriQuad({F X},W,G);
        else
            [TR,W,G,SMi] = QuadQuad({F X},W,G);
        end
        if i==1
            SM = SMi;
        else
            SM = SMi*SM;
        end
    end

    N1 = size(X,1);
    [F,X] = deal(TR{1},TR{2});
        
    % Project vertices onto unit sphere 
    N1 = N1 + 1;
    N2 = size(X,1);
    X(N1:N2,:) = ProjectOnSn(X(N1:N2,:));

end
TR = FormatOutputMesh(F,mean(R)*X,fmt);

