%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the averaged material parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% These arrays will store 1/epsr for x, y, and z-directed edges:
epsxm1 = ones(nx-1,ny,nz);
epsym1 = ones(nx,ny-1,nz);
epszm1 = ones(nx,ny,nz-1);


murxm1 = ones(nx,ny-1,nz-1);
murym1 = ones(nx-1,ny,nz-1);
murzm1 = ones(nx-1,ny-1,nz);


% X-average:
for i=1:nx-1
    for j=1:ny
        for k=1:nz
            epsr = 1;
            epsave = 0;
            numBlocks = 0;
            numBdries = 0;
            for b = 1:numOfDielBlocks
               if(i >= iEps1(b) && i < iEps2(b))
                   if(j >= jEps1(b) && j <= jEps2(b))
                       if(k >= kEps1(b) && k <= kEps2(b))
                           numBlocks = numBlocks+1;
                           epsave = epsave + epsrB(b);
                           
                           if(j == jEps1(b) || j == jEps2(b) ||...
                              k == kEps1(b) || k == kEps2(b) )
                              numBdries = numBdries+1;
                           end
                       end
                       %  Calculate 1/(average relative permittivity)
                       if(numBlocks == 0)
                           epsxm1(i,j,k) = 1.;
                       elseif(numBdries == 0 || numBdries == 4)
                           epsxm1(i,j,k) = 1.0/(epsave/numBlocks);
                       else
                           epsxm1(i,j,k) = 1.0/((epsave + 1.0)/numBlocks);
                       end
                   end
               end
            end
        end
    end
end


% Y-average:
for i=1:nx
    for j=1:ny-1
        for k=1:nz
            epsr = 1;
            epsave = 0;
            numBlocks = 0;
            numBdries = 0;
            for b = 1:numOfDielBlocks
               if(i >= iEps1(b) && i <= iEps2(b))
                   if(j >= jEps1(b) && j < jEps2(b))
                       if(k >= kEps1(b) && k <= kEps2(b))
                           numBlocks = numBlocks+1;
                           epsave = epsave + epsrB(b);
                           
                           if(i == iEps1(b) || i == iEps2(b) ||...
                              k == kEps1(b) || k == kEps2(b) )
                              numBdries = numBdries+1;
                           end
                       end
                       %  Calculate 1/(average relative permittivity)
                       if(numBlocks == 0)
                           epsym1(i,j,k) = 1.;
                       elseif(numBdries == 0 || numBdries == 4)
                           epsym1(i,j,k) = 1.0/(epsave/numBlocks);
                       else
                           epsym1(i,j,k) = 1.0/((epsave + 1.0)/numBlocks);
                       end
                   end
               end
            end
        end
    end
end


% Z-average:
for i=1:nx
    for j=1:ny
        for k=1:nz-1
            epsr = 1;
            epsave = 0;
            numBlocks = 0;
            numBdries = 0;
            for b = 1:numOfDielBlocks
               if(i >= iEps1(b) && i <= iEps2(b))
                   if(j >= jEps1(b) && j <= jEps2(b))
                       if(k >= kEps1(b) && k < kEps2(b))
                           numBlocks = numBlocks+1;
                           epsave = epsave + epsrB(b);
                           
                           if(i == iEps1(b) || i == iEps2(b) ||...
                              j == jEps1(b) || j == jEps2(b) )
                              numBdries = numBdries+1;
                           end
                       end
                       %  Calculate 1/(average relative permittivity)
                       if(numBlocks == 0)
                           epszm1(i,j,k) = 1.;
                       elseif(numBdries == 0 || numBdries == 4)
                           epszm1(i,j,k) = 1.0/(epsave/numBlocks);
                       else
                           epszm1(i,j,k) = 1.0/((epsave + 1.0)/numBlocks);
                       end
                   end
               end
            end
        end
    end
end


%  Process PEC blocks
for b = 1:numOfPEC
    % X-directed edges:
    for i = iPEC1(b):iPEC2(b)-1
        for j = jPEC1(b):jPEC2(b)
            for k = kPEC1(b):kPEC2(b)
                epsxm1(i,j,k) = 0.0;
            end
        end
    end
    % Y-directed edges:
    for i = iPEC1(b):iPEC2(b)
        for j = jPEC1(b):jPEC2(b)-1
            for k = kPEC1(b):kPEC2(b)
                epsym1(i,j,k) = 0.0;
            end
        end
    end
    % Z-directed edges:
    for i = iPEC1(b):iPEC2(b)
        for j = jPEC1(b):jPEC2(b)
            for k = kPEC1(b):kPEC2(b)-1
                epszm1(i,j,k) = 0.0;
            end
        end
    end
end


