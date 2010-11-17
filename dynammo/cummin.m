function [y] = cummin(x)
% cummulative min function
% y(i) = min(x(1:i))

m = x(1);
y = x;
for i = 1 : length(x)
    if (m < x(i))
        y(i) = m;
    else 
        m = x(i);
    end
end

