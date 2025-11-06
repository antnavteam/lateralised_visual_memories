function [Xfinal, Yfinal, idx_realclick] = fillpath (X, Y, step)

if nargin==2
    step = 0.05;
end
% X = X - X(1);
% Y = Y - Y(1);

Xfinal = [];Yfinal = [];
idx_realclick =[];


for i = 1 : (length (X) - 1);
    
x = X - X(i);
y = Y - Y(i);

[th, r] = cart2pol (x(i+1), y(i+1));

lfill = 0;
j = 1;
xfill = []; yfill = []; 
idx=[];

    while lfill < r

        [xfill(j), yfill(j)] = pol2cart (th, lfill);
        
        % fill the real click idx 1 = real click, 0 = fillpoint
        if j==1; idx(j)=1; % the first is a real point
        else idx(j)=0;
        end
        
        lfill = lfill+ step;
        j = j+1;
    end  

xfill = xfill + X(i);
yfill = yfill + Y(i);


Xfinal = [Xfinal, xfill];
Yfinal = [Yfinal, yfill];
idx_realclick = [idx_realclick, idx];
end


