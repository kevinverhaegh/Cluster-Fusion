function output=ShockShellModelGeneralEuler(rho, R0, Z, M, t, r0t, ReacCross)

%input: for each specie rho(r,specie), R0(specie), Z(specie), M(specie), 
%t(time), r0t(r), ReacCross(specie, specie). ReacCross is a square matrix
%indicating which species can fuse with each other, and which type of
%fusion reaction this is (DD=1, DT=2, PB=3). In the case of a single specie,
%rho becomes a vector and R0, Z, M bec, ReacCross become a single value.

%Determine and print number of species
sizevec=size(rho);
NrSpecies=min(sizevec);


disp('Number of species:')
disp(NrSpecies)

%initialise zero matrix
N=zeros(NrSpecies,1);

%fill matrix with the amount of ions per specie
for i=1:NrSpecies
N(i)=trapz(r0t, rho(:,i).*4.*pi.*r0t.^2);
end

%display the amount of ions per specie
disp('Amount of ions:')
disp(N)

%calculate the average density per specie
densavg=sum(Z.*N)./(4/3 * pi .* R0.^3); 

%calculate the characteristic expansion times per specie and display the
%characteristic time
t0 = sqrt(3.*M.*PhysConst.epsilon0./(densavg .* Z * PhysConst.e^2));

disp('CE time scales:')
disp(t0)

%rescale the radius onto the maximum radius and dimensionless
s0=r0t./R0(i);

%rescale the time for each ion specie onto the characteristic time (and
%dimensionless)
tau=zeros(NrSpecies,numel(t));
for i=1:NrSpecies
tau(i,:) = t./t0(i);
end

%Calculate for each specie the charge each spherical shell contains
%cumtrapz is a cumulative trapezoidal integral, which yields the number of
%the enclosed charge at each t=0 for each spherical shell for each ion
%specie. The diff() of that function delivers the difference of that
%function for each time steps (and the length of diff(A) is X-1 if the
%length of A ix X. Therefore, there is a zero added to Q, which means that
%the first spherical shell (r0t(1)) does not contain any charge (e.q.
%r0t(1)=0). The ' is there to transpose the vector in the right direction.
Q=zeros(NrSpecies,numel(r0t));
for i=1:NrSpecies
Q(i,:)=[0 diff(cumtrapz(r0t, 4.*pi.*Z(i).*rho(:,i).*r0t.^2))'];
%Normalise on the charge of the full cluster (seperately for each ion specie)
Q(i,:)=Q(i,:)./(4/3 .* pi .* R0(i).^3 .* densavg(i));
end

%initialise a zero vector which contains (index of shell, index of time,
%index of atom specie), which will be filled with the complete dynamics of
%the expansion
rvecN=zeros(numel(rho(:,1)), numel(t),NrSpecies);

%fill the first time step of the vector rvecN with the begin condition of
%the differential equation
for i=1:NrSpecies
rvecN(:,1,i)=r0t;
end
%introduce a velocity vector, similarly to rvecN
vvecN=0.*rvecN;

%loop over each new timestep
for i=2:numel(t)
    %loop over each particle specie
    for k=1:NrSpecies
        %loop over each spherical shell
        for j=2:numel(s0)
            %find the charge enclosed by this spherical shell, see findQ.m
            Qi=findQ(Q,rvecN,i,j,k);
            %2d order Euler method of solving the differential equation
            rvecN(j,i,k)=rvecN(j,i-1,k) + (t(i)-t(i-1))*vvecN(j,i-1,k) + (R0(k)^3)*((tau(k,i)-tau(k,i-1))^2) * (Qi/(rvecN(j,i-1,k).^2)) ;
            vvecN(j,i,k)=(rvecN(j,i,k)-rvecN(j,i-1,k))/(t(i)-t(i-1)) ; 
        end
    end
end


%determine the amount of particles per shell per species (similarly to
%determining Q.
NumPVir=zeros(numel(r0t)+1,NrSpecies);
for i=1:NrSpecies
NumPVirH=cumtrapz([0; r0t], 4*pi*[rho(1,i); rho(:,i)].*[0; r0t].^2);
NumPVir(:,i)=[NumPVirH(1); diff(NumPVirH)];
end

%Determine how many fusable species are present (e.q. in a DT cluster -
%that is 2 because DD and DT fusion reactions can occur, in a DD cluster
%that is 1 because only DD fusion reactions can occur)
ReacSpecies = sum(sum(ReacCross~=0));

%initialise the required zero vectors/data structures to hold the
%information required for each fusion reaction
totcross=zeros(ReacSpecies,1);
icrossvec=cell(ReacSpecies,1);
tcross=cell(ReacSpecies,1);
v1cross=cell(ReacSpecies,1);
v2cross=cell(ReacSpecies,1);
rcross=cell(ReacSpecies,1);
Erel=cell(ReacSpecies,1);
reac=cell(ReacSpecies,1);
vrel =cell(ReacSpecies,1);
FusyT = zeros(ReacSpecies,1);
fusind=0;
FusY=0;

%investigate for each fusionable specie (e.q. in a DT cluster DD is
%investigated once and DT is investigated once)
for k=1:numel(ReacCross(:,1))
    for l=1:numel(ReacCross(:,1))
        %if fusion reactions can occur between specie k and specie l
        if ReacCross(k,l)~=0
            %detect fusion reactions
            fusind=fusind+1;
            %initialise cross value
            cross=0;
            %search for collisions, a collision is an intersection of the various lines
            %of various radii between the required species k and l (in the case of DD k=l). 
            %This is required to count the number of intersections - after
            %which the required data structures can be loaded
            for i=1:numel(s0) %for each radii
                 for j=1:numel(s0) %for each radii
                 icross = intersection(rvecN(i,:,k), rvecN(j,:,l)); %see if intersection occurs
                    if icross~=0 %if intersection occurs
                        totcross(fusind)=totcross(fusind)+1; %note down number of intersections
                    end
                 end
            end

            %for create zero vectors in the data structures with the size
            %of intersections occuring
            icrossvec{fusind,1}=zeros(totcross(fusind),1);
            tcross{fusind,1}=0.*icrossvec{fusind,1};
            v1cross{fusind,1}=0.*icrossvec{fusind,1};
            v2cross{fusind,1}=0.*icrossvec{fusind,1};
            rcross{fusind,1}=zeros(totcross(fusind),2);
            reac{fusind,1}=zeros(totcross(fusind),1);
            Erel{fusind,1}=zeros(totcross(fusind),1);
            vrel{fusind,1}=zeros(totcross(fusind),1);

            %search for collisions - this time not for counting the number
            %of collisions but for filling in the required data in the zero
            %vectors
            for i=1:numel(s0) %for each radii
                for j=1:numel(s0) %for each radii
                icross = intersection(rvecN(i,:,k), rvecN(j,:,l)); %search for intersection
                    if icross~=0 %if intersection occurs
                            cross=cross+1; %add 1 to the counter of intersections
                            icrossvec{fusind,1}(cross)=round(icross); 
                            %note down the index of time when the intersection occurs (and round this off to an integer)
                            rcross{fusind,1}(cross, :)=[i, j]; 
                            %note down the pair of index of radii where the intersection occurs. e.q. [2,3] indicates that shell 2 of specie k intersects with shell 3 of specie l (in the case of DD k=l)
                            tcross{fusind,1}(cross)=findv(t,icross); %note down the exact time at which the intersection occurs
                            v1cross{fusind,1}(cross)=findv(vvecN(i,:,k),icross); %notes down the velocity of the i shock shell of specie k
                            v2cross{fusind,1}(cross)=findv(vvecN(j,:,l),icross); %notes down the velocity of the j shock shell of specie l
                     end
                 end
             end

            %calculate the center of mass energy of the collision of
            %the two shells in keV (1000 * e): (0.5 * M(l) * |v_i - v_j|^2 / (1000 * e)
                
            Erel{fusind,1} = 0.5 * M(l) * abs(v1cross{fusind,1} - v2cross{fusind,1}).^2 ./ (1000*PhysConst.e);
            
            %Checks which type of fusion reaction is applicable
            
            if ReacCross(k,l)==1 %DD reaction found
                load('DD_crossection.mat'); %load DD cross-section table (sigma(Epsilon))
                for i=1:numel(reac{fusind,1}) %look at all collisions
                    if Erel{fusind,1}(i)<20 %energy too low for fusion to occur
                        reac{fusind,1}(i)=0;
                    else %energy high enough for fusion to occur
                        [~,r]=min(abs(Erel{fusind,1}(i)-EnDDsp)); 
                        %Determine the cross-section (EnDDsp is a value
                        %with energies on which the cross-section is
                        %evaluated, which is loaded from the load file
                        reac{fusind,1}(i)=DDsp(r)*(10^-28)*abs(v1cross{fusind,1}(i)-v2cross{fusind,1}(i)); 
                        %the cross-section DDsp(r) in barns (10^-28 m^2) is multiplied
                        %times the relative velocity to determine the
                        %reactivity
                    end
                end
             else if ReacCross(k,l)==2
                            load('DT_crossection.mat');
                                for i=1:numel(reac{fusind,1})
                                    if Erel{fusind,1}(i)<7
                                        reac{fusind,1}(i)=0;
                                    else
                                        [~,r]=min(abs(Erel{fusind,1}(i)-EnDTsp));
                                        reac{fusind,1}(i)=DTsp(r)*(10^-28)*abs(v1cross{fusind,1}(i)-v2cross{fusind,1}(i));
                                    end
                                end
                         else if ReacCross(k,l)==3
                            load('PB_crossection.mat');
                                for i=1:numel(reac{fusind,1})
                                    if Erel{fusind,1}(i)<56
                                        reac{fusind,1}(i)=0;
                                    else
                                        [~,r]=min(abs(Erel{fusind,1}(i)-EnPBsp));
                                        reac{fusind,1}(i)=PBsp(r)*(10^-28)*abs(v1cross{fusind,1}(i)-v2cross{fusind,1}(i));
                                    end
                                end
                             end
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
                    
                    vrel{fusind,1}(i) = v1cross{fusind,1}(i) - v2cross{fusind,1}(i);
                    if vrel{fusind,1}(i)<0
                        if rcross{fusind,1}(i,1)~=numel(rvecN(:,1,k)) && rcross{fusind,1}(i,2)~=numel(rvecN(:,1,k))
                            tdur(i)=(rvecN(rcross{fusind,1}(i,1)+1, icrossvec{fusind,1}(i),k) - rvecN(rcross{fusind,1}(i,1), icrossvec{fusind,1}(i),k))/vrel{fusind,1}(i);
                            else if rcross{fusind,1}(i,1)~=1
                                tdur(i)=(rvecN(rcross{fusind,1}(i,1), icrossvec{fusind,1}(i),k) - rvecN(rcross{fusind,1}(i,1)-1, icrossvec{fusind,1}(i),k))/vrel{fusind,1}(i);
                                end
                        end
                    else if vrel{fusind,1}(i)>0
                        if rcross{fusind,1}(i,1)~=numel(rvecN(:,2,k)) && rcross{fusind,1}(i,2)~=numel(rvecN(:,1,l))
                            tdur(i) = (rvecN(rcross{fusind,1}(i,2)+1, icrossvec{fusind,1}(i),l) - rvecN(rcross{fusind,1}(i,2), icrossvec{fusind,1}(i),l))/vrel{fusind,1}(i);
                            else if rcross{fusind,1}(i,1)~=1
                                tdur(i) =  (rvecN(rcross{fusind,1}(i,2), icrossvec{fusind,1}(i),l) - rvecN(rcross{fusind,1}(i,2)-1, icrossvec{fusind,1}(i),l))/vrel{fusind,1}(i);
                                end
                        end
                        end
                    end
                    
                    FusY(i)=abs((V1(i)+V2(i))*n1(i)*n2(i)*(tdur(i))*reac{fusind,1}(i));
                end
                FusyT(fusind)=sum(FusY);
            end
        end
end


disp('Total Fusion Yield per cluster')
disp(sum(FusyT));
output.FusY=sum(FusY);


output.rvecN=rvecN;
output.Q=Q;
output.vvecN=vvecN;
output.NumPVir=NumPVir;
if cross==0
    output.cross=cross;
else
    output.cross=cross;
    output.rcross=rcross;
    output.tcross=tcross;
    output.v1cross=v1cross;
    output.v2cross=v2cross;
end
end
