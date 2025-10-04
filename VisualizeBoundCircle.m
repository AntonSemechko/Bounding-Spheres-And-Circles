function H = VisualizeBoundCircle(X,R,C,Xb,ha)
% Visualize a 2D point cloud and its minimum bounding circle (MBC).
%
% INPUT
%   - X     : M-by-2 list of point co-ordinates 
%   - R     : radius of the circle
%   - C     : 1-by-2 vector specifying the centroid of the circle
%   - Xb    : (optional) K-by-2 array of MBC support points; 2<=K<=3. 
%   - ha    : (optional) handle of the axis where bounding circle is to be
%             visualized.
%
% OUTPUT
%   - H     : 1-by-2 vector containing handles for the following objects: 
%               H(1)    : handle of the point cloud
%               H(2)    : handle for the circle
%               H(3)    : points supporting the MBC
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<2 || isempty(R) || isempty(C)
    [R,C] = ExactMinBoundCircle(X);
end

% Construct bounding circle
t = linspace(0,2*pi,1E3);
x = R*cos(t) + C(1);
y = R*sin(t) + C(2);

% Visualize point cloud
if nargin<4
    Xb = [];
elseif ~isnumeric(Xb) || ~ismatrix(Xb) || size(Xb,2)~=2 || size(Xb,1)>3
    error('Invalid format for 4th input argument (Xb)')
end

if nargin<5 || isempty(ha) || numel(ha)>1 || ~ishandle(ha) || ~strcmp(get(ha,'type'),'axes')
    figure('color','w')
else
    axes(ha)
end

H = zeros(1,2);
H(1) = plot(X(:,1),X(:,2),'.k','MarkerSize',15);
hold on

% Visualize bounding circle
H(2) = plot(x,y,'--k','LineWidth',0.5);
axis equal off tight 

if ~isempty(Xb)
    H(3) = plot(Xb(:,1),Xb(:,2),'.b','MarkerSize',30);
end

if nargout<1, clear H; end
