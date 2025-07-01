%Ez1Out(n) = Ez(3,153,2);
%Ez2Out(n) = Ez(3,103,2);
%Ez2bOut(n) = Ez(3,83,2);
%Ez3Out(n) = Ez(3,63,2);
%Ez4Out(n) = Ez(3,43,2);
%Ez5Out(n) = Ez(3,23,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Output Data Probes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputting the data probes per time step

%  Loop over the output quantities:
for iout = 1:numOutputQty
    
    switch OutputType(iout)
        case 1     % Ex
            sum = 0.0;
            numFlds = 0;
            for k = kOutput1(iout):kOutput2(iout)
                for j = jOutput1(iout):jOutput2(iout)
                    for i = iOutput1(iout):iOutput2(iout)-1
                        numFlds = numFlds+1;
                        sum = sum + ex(i,j,k)/dx;      
                    end
                end
            end
            if(numFlds > 0)
                OutputValue(n,iout) = sum/numFlds;  % can average over multiple edges (with caution)
            end
            
        case 2     % Ey
            sum = 0.0;
            numFlds = 0;
            for k = kOutput1(iout):kOutput2(iout)
                for j = jOutput1(iout):jOutput2(iout)-1
                    for i = iOutput1(iout):iOutput2(iout)
                        numFlds = numFlds+1;
                        sum = sum + ey(i,j,k)/dy;
                    end
                end
            end
            if(numFlds > 0)
                OutputValue(n,iout) = sum/numFlds;  % can average over multiple edges (with caution)
            end
            
        case 3     % Ez
            sum = 0.0;
            numFlds = 0;
            for k = kOutput1(iout):kOutput2(iout)-1
                for j = jOutput1(iout):jOutput2(iout)
                    for i = iOutput1(iout):iOutput2(iout)
                        numFlds = numFlds+1;
                        sum = sum + ez(i,j,k)/dz;
                    end
                end
            end
            if(numFlds > 0)
                OutputValue(n,iout) = sum/numFlds;  % can average over multiple edges (with caution)
            end
            
            
        case 4     % Hx
            sum = 0.0;
            numFlds = 0;
            for k = kOutput1(iout):kOutput2(iout)-1
                for j = jOutput1(iout):jOutput2(iout)-1
                    for i = iOutput1(iout):iOutput2(iout)
                        numFlds = numFlds+1;
                        sum = sum + hx(i,j,k)/(dx*eta0);
                    end
                end
            end
            if(numFlds > 0)
                OutputValue(n,iout) = sum/numFlds;  % can average over multiple edges (with caution)
            end
            
            
        case 5     % Hy
            sum = 0.0;
            numFlds = 0;
            for k = kOutput1(iout):kOutput2(iout)-1
                for j = jOutput1(iout):jOutput2(iout)
                    for i = iOutput1(iout):iOutput2(iout)-1
                        numFlds = numFlds+1;
                        sum = sum + hy(i,j,k)/(dy*eta0);
                    end
                end
            end
            if(numFlds > 0)
                OutputValue(n,iout) = sum/numFlds;  % can average over multiple edges (with caution)
            end
            
            
        case 6     % Hz
            sum = 0.0;
            numFlds = 0;
            for k = kOutput1(iout):kOutput2(iout)
                for j = jOutput1(iout):jOutput2(iout)-1
                    for i = iOutput1(iout):iOutput2(iout)-1
                        numFlds = numFlds+1;
                        sum = sum + hz(i,j,k)/(dz*eta0);
                    end
                end
            end
            if(numFlds > 0)
                OutputValue(n,iout) = sum/numFlds;  % can average over multiple edges (with caution)
            end
            
            
        case 7     % Voltage
            sum = 0.0;
            numparallel = 0;
            if OutputAxis(iout) == 1 % voltage integral along x
                for k = kOutput1(iout):kOutput2(iout)
                    for j = jOutput1(iout):jOutput2(iout)
                        numparallel = numparallel+1;
                        for i = iOutput1(iout):iOutput2(iout)-1
                            sum = sum - ex(i,j,k);   
                        end
                    end
                end
            elseif OutputAxis(iout) == 2 % voltage integral along y
                for k = kOutput1(iout):kOutput2(iout)
                    for i = iOutput1(iout):iOutput2(iout)
                        numparallel = numparallel+1;
                        for j = jOutput1(iout):jOutput2(iout)-1
                            sum = sum - ey(i,j,k); 
                        end
                    end
                end
            elseif OutputAxis(iout) == 3 % voltage integral along z
                for j = jOutput1(iout):jOutput2(iout)
                    for i = iOutput1(iout):iOutput2(iout)
                        numparallel = numparallel+1;
                        for k = kOutput1(iout):kOutput2(iout)-1
                            sum = sum - ez(i,j,k); 
                        end
                    end
                end
            end
            if(numparallel > 0)
                OutputValue(n,iout) = sum/numparallel;  % can average over parallel line integrals (with caution)
            end
                       
            
        case 8     %  Current
            sum = 0.0;
            numparallel = 0;
            if OutputAxis(iout) == 1 % current axis is || to x
                for i = iOutput1(iout):iOutput2(iout)
                    numparallel = numparallel+1;
                    kup = kOutput2(iout);
                    kdown = kOutput1(iout)-1;
                    for j = jOutput1(iout):jOutput2(iout)
                        sum = sum + (hy(i,j,kup)-hy(i,j,kdown)); 
                    end
                    jleft = jOutput1(iout)-1;
                    jright = jOutput2(iout);
                    for k = kOutput1(iout):kOutput2(iout)
                        sum = sum + (hz(i,jleft,k)-hz(i,jright,k));  
                    end
                end
            elseif OutputAxis(iout) == 2 % current axis is || to y
                for j = jOutput1(iout):jOutput2(iout)
                    numparallel = numparallel+1;
                    iup = iOutput2(iout);
                    idown = iOutput1(iout)-1;
                    for k = kOutput1(iout):kOutput2(iout)
                        sum = sum + (hz(iup,j,k)-hz(idown,j,k)); 
                    end
                    kleft = kOutput1(iout)-1;
                    kright = kOutput2(iout);
                    for i = iOutput1(iout):iOutput2(iout)
                        sum = sum + (hx(i,j,kleft)-hx(i,j,kright));  
                    end
                end
            elseif OutputAxis(iout) == 3 % current axis is || to z
                for k = kOutput1(iout):kOutput2(iout)
                    numparallel = numparallel+1;
                    jup = jOutput2(iout);
                    jdown = jOutput1(iout)-1;
                    for i = iOutput1(iout):iOutput2(iout)
                        sum = sum + (hx(i,jup,k)-hx(i,jdown,k));  
                    end
                    ileft = iOutput1(iout)-1;
                    iright = iOutput2(iout);
                    for j = jOutput1(iout):jOutput2(iout)
                        sum = sum + (hy(ileft,j,k)-hy(iright,j,k));
                    end
                end
            end
            if(numparallel > 0)
                OutputValue(n,iout) = -sum/numparallel/eta0;  % can average over parallel line integrals (with caution)
                                                         % negative for
                                                         % +ve axis 
            end  
    end
    
end