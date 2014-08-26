function M=movi(input)

species = numel(input.AMUZcombo(:,1));

hfig=(figure);
set(hfig,'Position',[0,0,1280,1280]);

%M=zeros(numel(input.t));

for j=1:100:numel(input.t)
    try
    for i=1:species
        if i==1 
            x=unique(input.x{1,1}(j,:), 'stable');
            y=unique(input.y{1,1}(j,:), 'stable');
            z=unique(input.z{1,1}(j,:), 'stable');
            S=unique(input.x{1,1}(j,:), 'stable');
            if numel(S)~=1
                S=10*ones(numel(unique(input.x{1,1}(j,:), 'stable')),1);
            end
            C=unique(input.x{1,1}(j,:), 'stable');
            if numel(C)~=1
                C=ones(numel(unique(input.x{1,1}(j,:), 'stable')),1);
            end
        else
        x=addvec(x, unique(input.x{i,1}(j,:), 'stable'), 1);
        y=addvec(y, unique(input.y{i,1}(j,:), 'stable'), 1);
        z=addvec(z, unique(input.z{i,1}(j,:), 'stable'), 1);
        Sh=unique(input.x{i,1}(j,:), 'stable');
        if numel(Sh)~=1
                Sh=50*ones(numel(unique(input.x{i,1}(j,:), 'stable')),1);
        end
        S=addvec(S, Sh,1);
        Ch=unique(input.x{i,1}(j,:), 'stable');
        if numel(Ch)~=1
                Ch=i*ones(numel(unique(input.x{i,1}(j,:), 'stable')),1);
        end
        C=addvec(C, Ch,1);
        end
    end
    [~,R]=min(abs(S-50));
    axiss=1.2*max(x(R:end));
    scatter3(x,y,z,S,C,'fill')
    legend
    title(num2str(input.t(j)));
    axis([-axiss, axiss, -axiss, axiss, -axiss, axiss]);
    M((j-1)/10 + 1)=getframe();
    catch
        disp('Couldnt display frame')
        disp('index')
        disp(j)
        perror('Quitting program')
    end
end