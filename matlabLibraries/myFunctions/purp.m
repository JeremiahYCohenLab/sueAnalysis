function c = purp(m)
%   Shades of cyan and red color map
%   purp(M) returns an M-by-3 matrix containing a "purple" colormap.

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

if rem(m,2)
    r = (0:0.5*(m+1)-1)'/max(0.5*(m+1)-1,1);
    c1 = [ones(0.5*(m+1),1) zeros(0.5*(m+1),1) r];
    c2 = [1-r r ones(0.5*(m+1),1)];
    c = [c1; c2(2:end,:)];   
else
    r = (0:0.5*m-1)'/max(0.5*m-1,1);
    c1 = [ones(0.5*m,1) zeros(0.5*m,1) r];
    r = (0:0.5*m)'/max(0.5*m,1);
    c2 = [1-r r ones(0.5*m+1,1)];
    c = [c1; c2(2:end,:)];
end
     
