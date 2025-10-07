function TR = FormatOutputMesh(Tri,V,fmt)
%
% INPUT
%   - Tri   : face array
%   - V     : vertex array
%   - fmt   : output format. See GetMeshData.m for the list of supported
%             formats.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


switch fmt
    case 1
        TR = triangulation(Tri,V);
    case 2
        TR = TriRep(Tri,V); %#ok<*DTRIREP>
    case 3
        TR = {Tri V};
    case 4
        TR = struct('faces',Tri,'vertices',V);
    otherwise
        TR = V;
end
