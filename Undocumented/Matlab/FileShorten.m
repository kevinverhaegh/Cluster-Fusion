function short=FileShorten(input,fract)
len=numel(input);
newlen = floor(len/fract);
short=cell(newlen,1);

for i=1:1:newlen
    newlen{i,1}=input{i*fract};
end