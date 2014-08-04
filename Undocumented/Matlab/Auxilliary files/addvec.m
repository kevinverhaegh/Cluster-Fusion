function output=addvec(vec1, vec2, nonzero)
try
output=[vec1 vec2];
catch
    try
        output=[vec1; vec2;];
    catch
        try
            output=[vec1; vec2];
        catch
            try
                output=[vec1 vec2;];
            catch
                error('Something wrong with the input vectors')
            end
        end
    end
end

if nonzero==1
    if numel(vec1)==1
        output=vec2;
    end
    if numel(vec2)==1
        output=vec1;
    end
end