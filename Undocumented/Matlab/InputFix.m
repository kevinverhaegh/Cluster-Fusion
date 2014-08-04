function fix = InputFix(input)
%disp('Fixing non monotonous time array')
for i=2:numel(input)
    if input{i,1}.param.time < input{i-1,1}.param.time
        break;
    end
end

if numel(input) == i
    %disp('No problems found')
    fix = input;
    return
end

fix = cell(i-1,1);

for j=1:i-1
    fix{j,1} = input{j,1};
end