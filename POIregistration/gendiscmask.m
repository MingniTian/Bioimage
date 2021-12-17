function [discmask] = gendiscmask(N)
%GENDISCMASK
%   [discmask] = gendiscmask(N) generates a N*N matrix. Elements within the
%   central region of the matrix, which is a disc with a radius of (N-1)/2,
%   have value 1 and all other elements have value 0.

%   Written in 2020 by Buyun Tian @ Tao Xu's lab, Institute of Biophysics,
%   Beijing, China.

r = (N - 1)/2;
centre = (N + 1)/2;
discmask = zeros(N);
for i = 1:N
    for j = 1:N
        discmask(i,j) = (norm([j i] - [centre centre]) <= r);
    end
end