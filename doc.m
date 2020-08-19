%% planes_intersection
% Function to compute the intersection between P1(M1,n1) and P2(M2,n2) planes of the 3D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2020.
%
%% Syntax
%
% [I, u, rc] = planes_intersection(n1, M1, n2, M2);
%
% [I, u, rc] = planes_intersection(n1, M1, n2, M2, verbose);
%
%% Description
%
% [I, u, rc] = planes_intersection(n1, M1, n2, M2) compute the
% intersection between the planes P1 and P2 defined by (M1,n1) and (M2,n2).
% The value of the return code rc provides the nature of the intersection
% (line, plane, or void). When the intersection is a line, its description
% is given by the point I and the director vector u.
%
% [I, u, rc] = planes_intersection(n1, M1, n2, M2, verbose) enables
% verbose mode when verbose is set to logical *true or numeric *1,
% and disables it when it is set to logical false or numeric 0.
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/73760-line-plane-intersection-3d?s_tid=prof_contriblnk line_plane_intersection> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73922-lines-intersection-3d-2d?s_tid=prof_contriblnk lines_intersection> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73490-point-to-plane-distance?s_tid=prof_contriblnk point_to_plane_distance> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73478-point-to-line-distance-3d-2d?s_tid=prof_contriblnk point_to_line_distance>
%
%% Principle
% Descartes plane equation & 3D line parametric equation
%
% - (1) : line director vector is the cross product of n1 & n2.
%
% - (2) : looking for a point belonging to the line (common between (M1,n1)
% and (M2,n2) planes.
%
% Solving the system :
%
% * a1+x + b1*y + c1*z + d1 = 0 (P1)
%
% * a2+x + b2*y + c2*z + d2 = 0 (P2)
%
% with Kramer technique. 
%
%% Inputs arguments
%
% - M1 : real row or column vector double, a point belonging to P1. numel(M1) = 3.
%
% - n1 : real row or column vector double, one P1 normal vector. numel(n1) = 3.
%
% - M2 : real row or column vector double, a point belonging to P2. numel(M2) = 3.
%
% - n2 : real row or column vector double, one P2 normal vector. numel(n2) = 3.
%
% - verbose : either logical *true/false or numeric *1/0.
%
%% Output arguments
%
% - I : real vector double, one point belonging to the line intersection. size(I) = size(n1).
%
% - u : real vector double, one line director vector. size(u) = size(n1).
%
% - rc : return code, numeric integer in the set {1,2,3}.
%
% * 0 = void / [] intersection
% * 1 = line intersection,
% * 2 = plane intersection, P2 = P1
%
%        rc return code is necessary to distinguish between cases where P1 and P2 intersection
%        is a line and where it is the plane itself (P1 and P2 coincide).
%
%% Example #1
% Intersection is a line
n1 = [1 1 1];
M1 = ones(1,3)/3;
n2 = [-1 1 0];
M2 = [1 1 0];
[I, u, rc] = planes_intersection(n1, M1, n2, M2) % expected : u = k*[1 1 -2], rc = 1

%% Example #2
% Intersection is void : P1 // P2
n1 = [0 0 1];
M1 = n1;
n2 = -n1;
M2 = n2;
[I, u, rc] = planes_intersection(n1, M1, n2, M2, true) % expected : I = [], u = [], rc = 0

%% Example #3
% Intersection is a plane : P2 = P1
n1 = [1 0 0];
M1 = n1;
n2 = [-2 0 0];
M2 = M1;
[I, u, rc] = planes_intersection(n1, M1, n2, M2, true) % expected : rc = 2