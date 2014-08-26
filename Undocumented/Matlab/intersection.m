function ind=intersection(s1, s2)
%finds the intersection between two vectors, assuming both vectors have the
%same size, and assuming that the vectors are not intersecting at the first
%or last index, it returns the fractional index of the intersection -
%depending on the linear interpolation near the intersection.

ind=0;
[~,r]=min(abs(s1(2:end-1) - s2(2:end-1))); %looks for the nearest point of the two vectors
r=r+1;
if (s1(r)==s2(r)) && (s1(r-1)~=s2(r-1)) && (s1(r+1) ~= s2(r+1))
    %the two vectors intersect exactly at point $r$, intersection found
    ind=r;
else if s1(r)==s2(r)
        %disp('lines are too similar, intersection occurs at multiple points')
        %Intersection occurs at multiple points - cannot return ind
    end
end

if ((s2(r)>s1(r)) && (s2(r+1)<s1(r+1))) || ((s2(r)<s1(r)) && (s2(r+1)>s1(r+1)))
    %intersection found - but the intersection point is not exactly at the
    %closest index, trying to find fractional index (e.q. ind = 2.25, would
    %imply that intersection occurs at 1/4th in between the second and
    %third point).
    %draw linear line through both trends
    a1=(s1(r+1) - s1(r));
    b1 = s1(r) - a1*(r);
    a2=(s2(r+1) - s2(r));
    b2 = s2(r) - a2*(r);
    ind=(b2-b1)/(a1-a2);

    else if (s2(r)>s1(r)) && (s2(r-1)<s1(r-1)) || ((s2(r)<s1(r)) && (s2(r-1)>s1(r-1)))
        %intersection found - but the intersection point is not exactly at the
        %closest index, trying to find fractional index (e.q. ind = 2.25, would
        %imply that intersection occurs at 1/4th in between the second and
        %third point).
        %draw linear line through both trends
        a1=(s1(r) - s1(r-1));
        b1 = s1(r) - a1*(r);
        a2=(s2(r) - s2(r-1));
        b2 = s2(r) - a2*(r);
        ind=(b2-b1)/(a1-a2);
        end
end
