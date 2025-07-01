%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inject Magnetic Field Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = n*dt;
coef_h = dt/mu0;

for is = 1:numSources
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the source signature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    srcA = 1.;
    if(sourceType(is) ~= 3)
        srcA = exp(-((t - srct_o(is))/srct_w(is))^2);
        if(sourceSigType(is) == 2)
            srcA = -((t - srct_o(is))/srct_w(is))*srcA;
        elseif(sourceSigType(is) == 3)
            srcA = sin(2*pi*srcfreq(is))*srcA;
        end
    end
    srcA = samp(is)*srcA;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % M-current source
    if(sourceType(is) ==2)
        % X-directed current source
        if(polInc(is) == 1)
            Mx_t = srcA*coef_h;
            for k = ksrc1(is):ksrc2(is)-1
                for j = jsrc1(is):jsrc2(is)-1
                    for i = isrc1(is):isrc2(is)
                        hx(i,j,k) = hx(i,j,k)-Mx_t*dx*eta0;
                    end
                end
            end
        elseif(polInc(is) == 2)
            My_t = srcA*coef_h;
            for k = ksrc1(is):ksrc2(is)-1
                for j = jsrc1(is):jsrc2(is)
                    for i = isrc1(is):isrc2(is)-1
                        hy(i,j,k) = hy(i,j,k)-My_t*dy*eta0;
                    end
                end
            end
        else
            Mz_t = srcA*coef_h;
            for k = ksrc1(is):ksrc2(is)
                for j = jsrc1(is):jsrc2(is)-1
                    for i = isrc1(is):isrc2(is)-1
                        hz(i,j,k) = hz(i,j,k)-Mz_t*dz*eta0;
                    end
                end
            end
        end
    end
    
end
