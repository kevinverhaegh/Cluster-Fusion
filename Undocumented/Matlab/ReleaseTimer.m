function ind = ReleaseTimer(input,points, thres)
n=0;
for i=1:numel(input)
    if abs(input(i))<thres
    n=n+1;
    else
        n=0;
    end
    if n==points
        ind=i-points+1;
        return
    end
end
ind = numel(input);