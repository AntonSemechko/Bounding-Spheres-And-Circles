function [TR, W, G, SM, idx_unq] = TriQuad(TR, W, G)
% Subdivide triangular surface mesh using generalized triangular 
% quadrisection. Triangular quadrisection is a linear subdivision procedure
% which inserts new vertices at the edge midpoints of the input mesh, 
% thereby producing four new faces for every face of the original mesh:
% 
%                     x3                        x3
%                    /  \      subdivision     /  \
%                   /    \        ====>       v3__v2
%                  /      \                  / \  / \
%                x1________x2              x1___v1___x2
%
%                      Original vertices : x1, x2, x3
%
%                      New vertices      : v1, v2, v3
%
%                      New faces         : [x1 v1 v3; x2 v2 v1; x3 v3 v2; v1 v2 v3] 
%
% In case of generalized triangular quadrisection, positions of the newly
% inserted vertices do not have to correspond to the edge midpoints, and 
% may be varied by assigning (positive) weights to the vertices of the 
% original mesh. For example, let xi and xj be two vertices connected by an
% edge, and suppose that Wi and Wj are the corresponding vertex weights. 
% Position of the new point on the edge (xi,xj) is defined as 
% (Wi*xi + Wj*xj)/(Wi+Wj). Note, ALL WEIGHTS MUST BE POSITIVE to avoid 
% degeneracies and self-intersections. 
%
% INPUT
%   - TR   : surface mesh represented as an object of 'TriRep' class,
%            object of 'triangulation' class, structure with 'faces' and
%            'vertices' fields, or a cell such that TR = {Tri,X}, where Tri
%            is a M-by-3 array of faces and X is a N-by-3 array of vertex
%            coordinates.
%   - W    : (optional) N-by-1 array of STRICTLY POSITIVE vertex weights 
%            used during interpolation of the new vertices, where N is the
%            total number of the original mesh vertices. 
%   - G    : (optional) scalar, vector, or tensor field defined at the mesh 
%            vertices and you wish to interpolate along with the vertices.
%            Note linear interpolation of tensor fields may produce output
%            that is not a valid tensor.
%
% OUTPUT
%   - TR  : subdivided mesh. Same format as the input mesh.
%   - W   : interpolated vertex weights.
%   - G   : interpolated scalar, vector, or tensor field defined at the
%           vertices of the subdivided mesh. Note, linear interpolation of
%           tensor fields may produce output that is not a valid tensor.
%   - SM  : sparse K-by-N subdivision matrix, where K is the number of
%           vertices in the subdivided mesh.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Get the list of vertex co-ordinates and list of faces
[Tri,X,fmt] = GetMeshData(TR);

% Make sure that the mesh is composed entirely of triangles 
if size(Tri,2)~=3
    error('This function is meant for triangular surface meshes, not quad or tet meshes')
end

% Check vertex weights
if nargin<2 || isempty(W)
    W = [];
elseif ~ismatrix(W) || numel(W)~=size(X,1) || sum(W<=eps)>0
    error('W must be a N-by-1 array with positive entries, where N is the # of mesh vertices')
else
    W = W(:);
end

% Field?
if nargin<3 || isempty(G)
    G = [];
elseif nargin==3 && ~isempty(G) && (size(G,1)~=size(X,1) || ~isnumeric(G) || ndims(G)>3)
    error('3rd input argument must be a %u-by-d array where d>=1',size(X,1))
end


% Edges
E = [Tri(:,1) Tri(:,2); Tri(:,2) Tri(:,3); Tri(:,3) Tri(:,1)];
E = sort(E,2);
[E,idx_unq,idx] = unique(E,'rows','stable'); % setOrder = 'stable' ensures that identical results will be obtained for meshes with the same connectivity

% Compute new vertex positions
if ~isempty(W) % insert new vertices based on vertex weights
    
    w = bsxfun(@rdivide,[W(E(:,1)),W(E(:,2))], W(E(:,1)) + W(E(:,2)));
    V = bsxfun(@times,X(E(:,1),:),w(:,1)) + bsxfun(@times,X(E(:,2),:),w(:,2));
    
    if ~isempty(G) && nargout>2
        G = cat(1, G, bsxfun(@times,G(E(:,1),:,:),w(:,1)) + bsxfun(@times,G(E(:,2),:,:),w(:,2)));
    end
    
    if nargout>1
        W = [W; W(E(:,1)).*w(:,1) + W(E(:,2)).*w(:,2)];
    end
        
else % insert new vertices at the edge mid-points
    
    V = (X(E(:,1),:) + X(E(:,2),:))/2;   
    if ~isempty(G) && nargout>2
        G = cat(1, G, (G(E(:,1),:,:) + G(E(:,2),:,:))/2);
    end
    
end

% Generate a subdivision matrix if one is required
Nx = size(X,1);   % # of vertices
Nf = size(Tri,1); % # of faces
if nargout>3
    
    Ne = size(E,1);
    if isempty(W), w = repmat([1 1]/2,[Ne 1]); end

    i = (1:Ne)' + Nx;
    i = cat(1,(1:Nx)',i,i);  

    j = cat(1,(1:Nx)',E(:));    
    w = cat(1,ones(Nx,1),w(:));
    
    SM = sparse(i, j, w, Nx + Ne, Nx, 2*Ne + Nx);    
    
end

% Assign indices to new triangle vertices
Vid1 = Nx + idx(1:Nf);
Vid2 = Nx + idx((Nf+1):2*Nf);
Vid3 = Nx + idx((2*Nf+1):3*Nf);

% Connectivities of the new faces
Tri1 = [Tri(:,1) Vid1  Vid3]; % [M x 3]
Tri2 = [Tri(:,2) Vid2  Vid1];
Tri3 = [Tri(:,3) Vid3  Vid2];
Tri4 = [Vid1     Vid2  Vid3];

Tri = reshape([Tri1 Tri2 Tri3 Tri4]',3,[])'; % [(4*3) x M] --> [(4*M) x 3]

% New mesh
X = [X; V];
TR = FormatOutputMesh(Tri,X,fmt);
