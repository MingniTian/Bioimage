function [cprocessed] = normxcorr2post(c, templatesize)
%NORMXCORR2POST
%   [cprocessed] = normxcorr2post(c, templatesize) post-processes c, the
%   output of normxcorr2post(template, A), based on templatesize such that 
%   the values computed with zero-padded edges are set to -1.

%   Written in 2020 by Buyun Tian @ Tao Xu's lab, Institute of Biophysics,
%   Beijing, China.

csize = size(c);
cprocessed = c;
cprocessed(1:(templatesize(1) - 1),:) = -1;
cprocessed((csize(1) - templatesize(1) + 2):csize(1),:) = -1;
cprocessed(:,1:(templatesize(2) - 1)) = -1;
cprocessed(:,(csize(2) - templatesize(2) + 2):csize(2)) = -1;
end