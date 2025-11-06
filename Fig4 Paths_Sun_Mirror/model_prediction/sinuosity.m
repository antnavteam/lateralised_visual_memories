function [Chunksvector_length mean_cumulated_angle bearings X_chunklocation, Y_chunklocation cumulated_angle] = sinuosity (X, Y, chunksize)
%
% X and Y are never modify, x and y are the chunk begning at 0.
bearings = [];

Pathend = 0;

refindex = 1; %index of the begning of the new chunk
endindex = length (X);% index of the end of the path
x = X; y = Y;
absoluteindex =1; %to keep track relative to the whole path

i = 1;
while Pathend == 0
    
    
    pathportion = [refindex:1:endindex]; %index of the path starting at the new chunk
    
    
    x = x (pathportion);
    y = y (pathportion);

    x = x - x(1); %set the begning of the new chunk at 0
    y = y - y(1);

    [th, r] = cart2pol (x,y);%each point in polar coordinate relative to the begining of th chunk
    
    a = find (r > chunksize); %all the point outside the chunksize circle

    if isempty (a);
        Pathend = 1;
    else
    
    bearings (i,1) = th (a(1)); %store the bearing of the last point of the current chunk
    
    refindex = a(1); %set the begining ref of the next chunk (=last point of current chunk)
    endindex = length (x); % set the end ref of the path left
    absoluteindex(i) = length (X)+1 - length (pathportion);%keep track of the chunk points
    i = i+1;
    end
end

X_chunklocation = X(absoluteindex)';
Y_chunklocation = Y(absoluteindex)';

if length(bearings)>1
Chunksvector_length = circ_r (bearings); % resulant length of the vector from each bearing (0 = loop; 1 = straight)
    
for j = 1 : (length (bearings) - 1);
    angle = sqrt ((bearings(j) - bearings (j+1))^2);
    if angle > pi;
        angle = angle - pi;
    end

   cumulated_angle(j) = angle;
end
mean_cumulated_angle = mean (cumulated_angle);


else
mean_cumulated_angle=NaN;
Chunksvector_length=NaN;
bearings=NaN;
end

