function out=evalu(y, x, x0)
[~,i] = min(abs(x - x0));
if x(i) - x0 > 0
    if i==1
        out=0;
        return
    end
    x1 = x(i - 1);
    x2 = x(i);
    y1 = y(i - 1);
    y2 = y(i);
end
if x(i) - x0 < 0
        x1 = x(i);
        if i==numel(x)
            out=0;
            return
        end
        x2 = x(i+1);
        y1 = y(i);
        y2 = y(i+1);
end
if x(i) == x0
    out = y(i);
    return
end

a = (y2 - y1)/(x2 - x1);
b = y2 - a*x2;

out = a*x0 + b;