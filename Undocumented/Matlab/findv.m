function val=findv(vec,ind)
%auxiliary function used to find the value (val) of an (fractional) index (ind)
%in a vector (vec)

%if the fractional index is an integer
if isinteger(ind)==1
    val=vec(ind);%return value of that index
else
     %if the fractional index is not an integer, perform linear interpolation on
     %the two values this index is in between. E.q. if ind = 2.25, draw a
     %linear line between vec(2) and vec(3) in order to obtain the value of
     %vec at 2.25. 
     lowind=floor(ind); %round off to lowest full integer (e.q. floor(2.999) = 2)
     a=(vec(lowind+1) - vec(lowind)); %determine slope of the linear interpolation
     b = vec(lowind) - a*(lowind); %determine offset of the linear interpolation
     val=a*ind + b; %determine value
end