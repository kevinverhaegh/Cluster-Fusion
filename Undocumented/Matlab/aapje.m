for i=1:numel(reac)
    if rcross(i,1)~=numel(rvecN(:,1))
    V1(i)=abs(4*pi*(rvecN(rcross(i,1), icrossvec(i)).^2) * (rvecN(rcross(i,1)+1, icrossvec(i)) - rvecN(rcross(i,1), icrossvec(i))));
    else if rcross(i,1)~=1
        V1(i)=abs(4*pi*(rvecN(rcross(i,1), icrossvec(i)).^2) * (rvecN(rcross(i,1), icrossvec(i)) - rvecN(rcross(i,1)-1, icrossvec(i))));
        end
    end
    
    if rcross(i,2)~=numel(rvecN(:,2))
    V2(i)=abs(4*pi*(rvecN(rcross(i,2), icrossvec(i)).^2) * (rvecN(rcross(i,2)+1, icrossvec(i)) - rvecN(rcross(i,2), icrossvec(i))));
    else if rcross(i,1)~=1
        V2(i)=abs(4*pi*(rvecN(rcross(i,2), icrossvec(i)).^2) * (rvecN(rcross(i,2), icrossvec(i)) - rvecN(rcross(i,2)-1, icrossvec(i))));
        end
    end
    
    n1(i)=NumPVir(rcross(i,1))/V1(i);
    n2(i)=NumPVir(rcross(i,2))/V2(i);
    if i~=numel(reac)
        if rcross(i,1)==rcross(i+1,1)
            tdur(i)=tcross(i+1)-tcross(i);
        else
            tdur(i)=0;
        end
    else
        tdur(i)=0;
    end
    FusY(i)=(V1(i)+V2(i))*n1(i)*n2(i)*tdur(i)*reac(i);
end
