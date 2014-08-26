ReacSpecies = sum(sum(ReacCross~=0));
totcross=zeros(ReacSpecies,1);
icrossvec=cell(ReacSpecies,1);
tcross=cell(ReacSpecies,1);
v1cross=cell(ReacSpecies,1);
v2cross=cell(ReacSpecies,1);
rcross=cell(ReacSpecies,1);
Erel=cell(ReacSpecies,1);
reac=cell(ReacSpecies,1);
FusyT = zeros(ReacSpecies,1);
fusind=0;
FusY=0;
for k=1:numel(ReacCross(:,1))
    for l=1:numel(ReacCross(:,1))
        if ReacCross(k,l)~=0
            %detect fusion reactions
            fusind=fusind+1;
            if k==l
                %same specie fusion
            else
                %different specie fusion
                cross=0;
                %search for collisions, a collision is an intersection of the various lines
                %of various radia. 
                for i=1:numel(s0)
                    for j=1:numel(s0)
                        icross = intersection(rvecN(i,:,k), rvecN(j,:,l));
                        if icross~=0
                            totcross(fusind)=totcross(fusind)+1;
                        end
                    end
                end

                icrossvec{fusind,1}=zeros(totcross,1);
                tcross{fusind,1}=0.*icrossvec{fusind,1};
                v1cross{fusind,1}=0.*icrossvec{fusind,1};
                v2cross{fusind,1}=0.*icrossvec{fusind,1};
                rcross{fusind,1}=zeros(totcross,2);
                reac{fusind,1}=zeros(totcross,1);
                Erel{fusind,1}=zeros(totcross,1);

                for i=1:numel(s0)
                    for j=1:numel(s0)
                        icross = intersection(rvecN(i,:,k), rvecN(j,:,l));
                        if icross~=0
                            cross=cross+1;
                            icrossvec{fusind,1}(cross)=round(icross);
                            rcross{fusind,1}(cross, :)=[i, j];
                            tcross{fusind,1}(cross)=findv(t,icross);
                            v1cross{fusind,1}(cross)=findv(vvecN(i,:,k),icross);
                            v2cross{fusind,1}(cross)=findv(vvecN(j,:,l),icross);
                        end
                    end
                end

                %calculate fusion reactions
                %load DD crossections
                    if ReacCross(k,l)==1
                    load('DD_crossection.mat');
                    end
                Erel{fusind,1} = 0.5 * M(l) * abs(v1cross{fusind,1} - v2cross{fusind,1}).^2 ./ PhysConst.e;

                for i=1:numel(reac{fusind,1})
                    if Erel{fusind,1}(i)<20
                        reac{fusind,1}(i)=0;
                    else
                        [~,r]=min(abs(Erel{fusind,1}(i)-EnDDsp));
                        reac{fusind,1}(i)=DDsp(r)*(10^-28)*abs(v1cross{fusind,1}(i)-v2cross{fusind,1}(i));
                    end
                end

                V1=reac{fusind,1}.*0;
                V2=reac{fusind,1}.*0;
                n1=0.*V1;
                n2=0.*V2;
                tdur=0.*V1;
                FusY=0.*V1;
                
                for i=1:numel(reac{fusind,1})
                    if rcross{fusind,1}(i,1)~=numel(rvecN(:,1,k)) && rcross{fusind,1}(i,2)~=numel(rvecN(:,1,k))
                    V1(i)=abs(4*pi*(rvecN(rcross{fusind,1}(i,1), icrossvec{fusind,1}(i),k).^2) * (rvecN(rcross{fusind,1}(i,1)+1, icrossvec{fusind,1}(i),k) - rvecN(rcross{fusind,1}(i,1), icrossvec{fusind,1}(i),k)));
                    else if rcross{fusind,1}(i,1)~=1
                        V1(i)=abs(4*pi*(rvecN(rcross{fusind,1}(i,1), icrossvec{fusind,1}(i),k).^2) * (rvecN(rcross{fusind,1}(i,1), icrossvec{fusind,1}(i),k) - rvecN(rcross{fusind,1}(i,1)-1, icrossvec{fusind,1}(i),k)));
                        end
                    end
    
                    if rcross{fusind,1}(i,1)~=numel(rvecN(:,2,k)) && rcross{fusind,1}(i,2)~=numel(rvecN(:,1,l))
                        V2(i)=abs(4*pi*(rvecN(rcross{fusind,1}(i,2), icrossvec{fusind,1}(i),l).^2) * (rvecN(rcross{fusind,1}(i,2)+1, icrossvec{fusind,1}(i),l) - rvecN(rcross{fusind,1}(i,2), icrossvec{fusind,1}(i),l)));
                    else if rcross{fusind,1}(i,1)~=1
                        V2(i)=abs(4*pi*(rvecN(rcross{fusind,1}(i,2), icrossvec{fusind,1}(i),l).^2) * (rvecN(rcross{fusind,1}(i,2), icrossvec{fusind,1}(i),l) - rvecN(rcross{fusind,1}(i,2)-1, icrossvec{fusind,1}(i),l)));
                        end
                    end
    
                    n1(i)=NumPVir(rcross{fusind,1}(i,1))/V1(i);
                    n2(i)=NumPVir(rcross{fusind,1}(i,2))/V2(i);
                    if i~=numel(reac{fusind,1})
                        if rcross{fusind,1}(i,1)==rcross{fusind,1}(i+1,1)
                            tdur(i)=tcross{fusind,1}(i+1)-tcross{fusind,1}(i);
                        else
                            tdur(i)=0;
                        end
                    else
                        tdur(i)=0;
                    end
                    FusY(i)=(V1(i)+V2(i))*n1(i)*n2(i)*tdur(i)*reac{fusind,1}(i);
                end
                FusyT(fusind)=sum(FusY);
            end
        end
    end
end

disp('Total Fusion Yield per cluster')
disp(sum(FusyT));
output.FusY=sum(FusY);