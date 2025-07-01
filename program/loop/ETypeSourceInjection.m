%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inject Electric Field Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (n-0.5)*dt;
coef_e = codt;

for is = 1:numSources
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the source signature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    srcA = 1.;
    if(sourceType(is) ~= 3)
        srcA = exp(-((t - srct_o(is))/srct_w(is))^2);
        if(sourceSigType(is) == 2)
            srcA = -(2.0*(t - srct_o(is))/srct_w(is))*srcA;
        elseif(sourceSigType(is) == 3)
            srcA = sin(2*pi*srcfreq(is))*srcA;
        end
    end
    srcA = samp(is)*srcA;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J-current source
    if(sourceType(is) == 1)
        % X-directed current source
        if(polInc(is) == 1)
            Jx_t = srcA*coef_e*eta0;
            for k = ksrc1(is):ksrc2(is)
                for j = jsrc1(is):jsrc2(is)
                    for i = isrc1(is):isrc2(is)-1
                        ex(i,j,k) = ex(i,j,k)-Jx_t*epsxm1(i,j,k)*dx;
                    end
                end
            end
        elseif(polInc(is) == 2)
            Jy_t = srcA*coef_e*eta0;
            for k = ksrc1(is):ksrc2(is)
                for j = jsrc1(is):jsrc2(is)-1
                    for i = isrc1(is):isrc2(is)
                        ey(i,j,k) = ey(i,j,k)-Jy_t*epsym1(i,j,k)*dy;
                    end
                end
            end
        else
            Jz_t = srcA*coef_e*eta0;
            ez(isrc1(is):isrc2(is),jsrc1(is):jsrc2(is),ksrc1(is):ksrc2(is)-1) = ez(isrc1(is):isrc2(is),jsrc1(is):jsrc2(is),ksrc1(is):ksrc2(is)-1)-Jz_t*epszm1(isrc1(is):isrc2(is),jsrc1(is):jsrc2(is),ksrc1(is):ksrc2(is)-1)*dz;
        end
        
    elseif(sourceType(is) == 3)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  E-PlaneWaveInjection  (TBD)
        error('Error:  plane wave source is not yet implemented');
        
    elseif(sourceType(is) == 4)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Delta-Gap Source  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Soft-voltage Source  (TBD)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Soft-voltage Source  
        % X-directed Soft-voltage Source
        switch polInc(is)
            case 1
                Vx_t = srcA/((isrc2(is)-isrc1(is)));
                ex(isrc1(is):isrc2(is)-1,jsrc1(is):jsrc2(is),ksrc1(is):ksrc2(is)) = -Vx_t;
        % Y-directed Soft-voltage Source
            case 2
                Vy_t = srcA/((jsrc2(is)-jsrc1(is)));
                ey(isrc1(is):isrc2(is),jsrc1(is):jsrc2(is)-1, ksrc1(is):ksrc2(is)) = -Vy_t;

            otherwise
        % Z-directed Soft-voltage Source
            Vz_t = srcA/((ksrc2(is)-ksrc1(is)));
            ez(isrc1(is):isrc2(is),jsrc1(is):jsrc2(is),ksrc1(is):ksrc2(is)-1) = -Vz_t;

        end
        
        
    elseif(sourceType(is) == 5)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Hard-voltage Source  (TBD)
        error('Error:  Hard-voltage source is not yet implemented');
        
    elseif(sourceType(is) == 6)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Soft-voltage Source  (TBD)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Soft-voltage Source  
        % X-directed Soft-voltage Source
        if(polInc(is) == 1)
            Vx_t = srcA/((isrc2(is)-isrc1(is)));
            for k = ksrc1(is):ksrc2(is)
                for j = jsrc1(is):jsrc2(is)
                    for i = isrc1(is):isrc2(is)-1
                        ex(i,j,k) = ex(i,j,k)-Vx_t;
                    end
                end
            end
        % Y-directed Soft-voltage Source
        elseif(polInc(is) == 2)
            Vy_t = srcA/((jsrc2(is)-jsrc1(is)));
            for k = ksrc1(is):ksrc2(is)
                for j = jsrc1(is):jsrc2(is)-1
                    for i = isrc1(is):isrc2(is)
                        ey(i,j,k) = ey(i,j,k)-Vy_t;
                    end
                end
            end
        else
        % Z-directed Soft-voltage Source
            Vz_t = srcA/((ksrc2(is)-ksrc1(is)));
            for k = ksrc1(is):ksrc2(is)-1
                for j = jsrc1(is):jsrc2(is)
                    for i = isrc1(is):isrc2(is)
                        ez(i,j,k) = ez(i,j,k)-Vz_t;
                    end
                end
            end
        end
        
    elseif(sourceType(is) == 7)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Thevenin Source  (TBD)
        error('Error:  Thevenin source is not yet implemented');
    end
    
end
