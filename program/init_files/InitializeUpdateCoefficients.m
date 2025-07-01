%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the update coefficients for e,h, qe and qh:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codt = c0*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e-update coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize:
coefex = ones(nx-1,ny,nz)*codt*dx/(dy*dz);
coefey = ones(nx,ny-1,nz)*codt*dy/(dz*dx);
coefez = ones(nx,ny,nz-1)*codt*dz/(dx*dy);


% scale by 1/epsr:
coefex = coefex.*epsxm1;
coefey = coefey.*epsym1;
coefez = coefez.*epszm1;

%%%%%%%%%%%%%%%%%%%%%%%%
% scale by kappa values:
%%%%%%%%%%%%%%%%%%%%%%%%
% Scale coefex by kappaX_h, b/c this is sampled at i+1/2, and
% normalize by kappaY_e and kappaZ_e, b/c these are sampled at j and k
for k = 2:nz-1
    for j = 2:ny-1
        for i = 1:nx-1
            coefex(i,j,k) = coefex(i,j,k)*kappaX_h(i)/(kappaY_e(j)*kappaZ_e(k));
        end
    end
end

% Scale coefey by kappaY_h, b/c this is sampled at j+1/2, and
% normalize by kappaZ_e and kappaX_e, b/c these are sampled at k and i
for k = 2:nz-1
    for j = 1:ny-1
        for i = 2:nx-1
            coefey(i,j,k) = coefey(i,j,k)*kappaY_h(j)/(kappaZ_e(k)*kappaX_e(i));
        end
    end
end

% Scale coefez by kappaZ_h, b/c this is sampled at k+1/2, and
% normalize by kappaX_e and kappaY_e, b/c these are sampled at i and j
for k = 1:nz-1
    for j = 2:ny-1
        for i = 2:nx-1
            coefez(i,j,k) = coefez(i,j,k)*kappaZ_h(k)/(kappaX_e(i)*kappaY_e(j));
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h-update coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize:
coefhx = ones(nx,ny-1,nz-1)*codt*dx/(dy*dz);
coefhy = ones(nx-1,ny,nz-1)*codt*dy/(dz*dx);
coefhz = ones(nx-1,ny-1,nz)*codt*dz/(dx*dy);

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  Thin wire coefficient modifications
% %%  H-field (assuming a z-directed wire)
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% for iw = 1:numThinWire
%  student fill in
% end


% scale by 1/mur:
%coefhx = coefhx.*murxm1;
%coefhy = coefhy.*murym1;
%coefhz = coefhz.*murzm1;

%%%%%%%%%%%%%%%%%%%%%%%%
% scale by kappa values:
%%%%%%%%%%%%%%%%%%%%%%%%
% Scale coefhx by kappaX_e, b/c this is sampled at i, and
% normalize by kappaY_h and kappaZ_h, b/c these are sampled at j+1/2 and k+1/2
for k = 1:nz-1
    for j = 1:ny-1
        for i = 1:nx
            coefhx(i,j,k) = coefhx(i,j,k)*kappaX_e(i)/(kappaY_h(j)*kappaZ_h(k));
        end
    end
end

% Scale coefhy by kappaY_e, b/c this is sampled at j, and
% normalize by kappaZ_h and kappaX_h, b/c these are sampled at k+1/2 and i+1/2
for k = 1:nz-1
    for j = 1:ny
        for i = 1:nx-1
            coefhy(i,j,k) = coefhy(i,j,k)*kappaY_e(j)/(kappaZ_h(k)*kappaX_h(i));
        end
    end
end

% Scale coefhz by kappaZ_e, b/c this is sampled at k, and
% normalize by kappaX_h and kappaY_h, b/c these are sampled at i+1/2 and j+1/2
for k = 1:nz
    for j = 1:ny-1
        for i = 1:nx-1
            coefhz(i,j,k) = coefhz(i,j,k)*kappaZ_e(k)/(kappaX_h(i)*kappaY_h(j));
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q-update coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PML Auxiliary Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%  PML E-Field-updates
%%%%%%%%%%%%%%%%
if(npmlx_1 > 0)
    b_eyz_xz_1 = zeros(npmlx_1+1,ny-1,nz-1);
    c_ey_xz_1 = zeros(npmlx_1+1,ny-1,nz-1);
    c_ez_xy_1 = zeros(npmlx_1+1,ny-1,nz-1);

    for k = 1:nz-1
        for j = 1:ny-1
            for i = 1:npmlx_1+1
                b_eyz_xz_1(i,j,k) = bx_e_1(i);
                c_ey_xz_1(i,j,k) = cx_e_1(i)*kappaY_h(j)*dy*codt...
                                   *epsym1(i,j,k)/(kappaZ_e(k)*dz);
                c_ez_xy_1(i,j,k) = cx_e_1(i)*kappaZ_h(k)*dz*codt...
                                   *epszm1(i,j,k)/(kappaY_e(j)*dy);
            end
        end
    end
end

if(npmlx_2 > 0)
    b_eyz_xz_2 = zeros(npmlx_2,ny-1,nz-1);
    c_ey_xz_2 = zeros(npmlx_2,ny-1,nz-1);
    c_ez_xy_2 = zeros(npmlx_2,ny-1,nz-1);
    
    for k = 1:nz-1
       for j = 1:ny-1
          i = nx - npmlx_2;
          for i2 = 1:npmlx_2
              
              b_eyz_xz_2(i2,j,k) = bx_e_2(i2);
              c_ey_xz_2(i2,j,k) = cx_e_2(i2)*kappaY_h(j)*dy*codt...
                                  *epsym1(i,j,k)/(kappaZ_e(k)*dz);
              c_ez_xy_2(i2,j,k) = cx_e_2(i2)*kappaZ_h(k)*dz*codt...
                                  *epszm1(i,j,k)/(kappaY_e(j)*dy);
             i = i + 1;
          end 
       end 
    end
end


if(npmly_1 > 0)
    b_ezx_yx_1 = zeros(nx-1,npmly_1+1,nz-1);
    c_ez_yx_1 = zeros(nx-1,npmly_1+1,nz-1);
    c_ex_yz_1 = zeros(nx-1,npmly_1+1,nz-1);

    for k = 1:nz-1
        for j = 2:npmly_1+1
            for i = 1:nx-1
                b_ezx_yx_1(i,j,k) = by_e_1(j);
                c_ez_yx_1(i,j,k) = cy_e_1(j)*kappaZ_h(k)*dz*codt...
                                   *epszm1(i,j,k)/(kappaX_e(i)*dx);
                c_ex_yz_1(i,j,k) = cy_e_1(j)*kappaX_h(i)*dx*codt...
                                   *epsxm1(i,j,k)/(kappaZ_e(k)*dz);
            end
        end
    end
end

if(npmly_2 > 0)
    b_ezx_yx_2 = zeros(nx-1,npmly_2,nz-1);
    c_ez_yx_2 = zeros(nx-1,npmly_2,nz-1);
    c_ex_yz_2 = zeros(nx-1,npmly_2,nz-1);
    
    for k = 1:nz-1
       j = ny - npmly_2;
       for j2 = 1:npmly_2
          for i = 1:nx-1
                b_ezx_yx_2(i,j2,k) = by_e_2(j2);
                c_ez_yx_2(i,j2,k) = cy_e_2(j2)*kappaZ_h(k)*dz*codt...
                                    *epszm1(i,j,k)/(kappaX_e(i)*dx);
                c_ex_yz_2(i,j2,k) = cy_e_2(j2)*kappaX_h(i)*dx*codt...
                                    *epsxm1(i,j,k)/(kappaZ_e(k)*dz);
          end
          j = j + 1;
       end 
    end
end

if(npmlz_1 > 0)
    b_exy_zy_1 = zeros(nx-1,ny-1,npmlz_1+1);
    c_ex_zy_1 = zeros(nx-1,ny-1,npmlz_1+1);
    c_ey_zx_1 = zeros(nx-1,ny-1,npmlz_1+1);

    for k = 2:npmlz_1+1
        for j = 1:ny-1
            for i = 1:nx-1
                b_exy_zy_1(i,j,k) = bz_e_1(k);
                c_ex_zy_1(i,j,k) = cz_e_1(k)*kappaX_h(i)*dx*codt...
                                   *epsxm1(i,j,k)/(kappaY_e(j)*dy);
                c_ey_zx_1(i,j,k) = cz_e_1(k)*kappaY_h(j)*dy*codt...
                                   *epsym1(i,j,k)/(kappaX_e(i)*dx);
            end
        end
    end
end

if(npmlz_2 > 0)
    b_exy_zy_2 = zeros(nx-1,ny-1,npmlz_1);
    c_ex_zy_2 = zeros(nx-1,ny-1,npmlz_1);
    c_ey_zx_2 = zeros(nx-1,ny-1,npmlz_1);
    
    k = nz - npmlz_2;
    for k2 = 1:npmlz_2
       for j = 1:ny-1
          for i = 1:nx-1
                b_exy_zy_2(i,j,k2) = bz_e_2(k2);
                c_ex_zy_2(i,j,k2) = cz_e_2(k2)*kappaX_h(i)*dx*codt...
                                    *epsxm1(i,j,k)/(kappaY_e(j)*dy);
                c_ey_zx_2(i,j,k2) = cz_e_2(k2)*kappaY_h(j)*dy*codt...
                                    *epsym1(i,j,k)/(kappaX_e(i)*dx);
          end 
       end 
       k = k + 1;
    end
end


%%%%%%%%%%%%%%%%
%  PML H-Field-updates
%%%%%%%%%%%%%%%%

if(npmlx_1 > 0)
    b_hyz_xz_1 = zeros(npmlx_1,ny-1,nz-1);
    c_hy_xz_1 = zeros(npmlx_1,ny-1,nz-1);
    c_hz_xy_1 = zeros(npmlx_1,ny-1,nz-1);
    for k = 1:nz-1
       for j = 1:ny-1
          for i = 1:npmlx_1
             b_hyz_xz_1(i,j,k) = bx_h_1(i);
             c_hy_xz_1(i,j,k) = cx_h_1(i)*kappaY_e(j)*dy*codt...
                                *murym1(i,j,k)/(kappaZ_h(k)*dz);
             c_hz_xy_1(i,j,k) = cx_h_1(i)*kappaZ_e(k)*dz*codt...
                                *murzm1(i,j,k)/(kappaY_h(j)*dy);
          end
       end
    end
end

if(npmlx_2 > 0)
    b_hyz_xz_2 = zeros(npmlx_2,ny-1,nz-1);
    c_hy_xz_2 = zeros(npmlx_2,ny-1,nz-1);
    c_hz_xy_2 = zeros(npmlx_2,ny-1,nz-1);
    
for k = 1:nz-1
   for j = 1:ny-1
      i = nx - npmlx_2;
      for i2 = 1:npmlx_2
             b_hyz_xz_2(i2,j,k) = bx_h_2(i2);
             c_hy_xz_2(i2,j,k) = cx_h_2(i2)*kappaY_e(j)*dy*codt...
                                *murym1(i,j,k)/(kappaZ_h(k)*dz);
             c_hz_xy_2(i2,j,k) = cx_h_2(i2)*kappaZ_e(k)*dz*codt...
                                *murzm1(i,j,k)/(kappaY_h(j)*dy);
         i = i + 1;
      end
   end
end
end


if(npmly_1 > 0)
    b_hzx_yx_1 = zeros(nx-1,npmly_1,nz-1);
    c_hz_yx_1 = zeros(nx-1,npmly_1,nz-1);
    c_hx_yz_1 = zeros(nx-1,npmly_1,nz-1);
    for k = 1:nz-1
       for j = 1:npmly_1
          for i = 1:nx-1
             b_hzx_yx_1(i,j,k) = by_h_1(j);
             c_hz_yx_1(i,j,k) = cy_h_1(j)*kappaZ_e(k)*dz*codt...
                                *murzm1(i,j,k)/(kappaX_h(i)*dx);
             c_hx_yz_1(i,j,k) = cy_h_1(j)*kappaX_e(i)*dx*codt...
                                *murxm1(i,j,k)/(kappaZ_h(k)*dz);
          end
       end
    end
end

if(npmly_2 > 0)
    b_hzx_yx_2 = zeros(nx-1,npmly_2,nz-1);
    c_hz_yx_2 = zeros(nx-1,npmly_2,nz-1);
    c_hx_yz_2 = zeros(nx-1,npmly_2,nz-1);
    for k = 1:nz-1
       j = ny - npmly_2;
       for j2 = 1:npmly_2
          for i = 1:nx-1
             b_hzx_yx_2(i,j2,k) = by_h_2(j2);
             c_hz_yx_2(i,j2,k) = cy_h_2(j2)*kappaZ_e(k)*dz*codt...
                                *murzm1(i,j,k)/(kappaX_h(i)*dx);
             c_hx_yz_2(i,j2,k) = cy_h_2(j2)*kappaX_e(i)*dx*codt...
                                *murxm1(i,j,k)/(kappaZ_h(k)*dz);
          end
          j = j + 1;
       end
    end
end

if(npmlz_1 > 0)
    b_hxy_zy_1 = zeros(nx-1,ny-1,npmlz_1);
    c_hx_zy_1 = zeros(nx-1,ny-1,npmlz_1);
    c_hy_zx_1 = zeros(nx-1,ny-1,npmlz_1);
    for k = 1:npmlz_1
       for j = 1:ny-1
          for i = 1:nx-1
             b_hxy_zy_1(i,j,k) = bz_h_1(k);
             c_hx_zy_1(i,j,k) = cz_h_1(k)*kappaX_e(i)*dx*codt...
                                *murxm1(i,j,k)/(kappaY_h(j)*dy);
             c_hy_zx_1(i,j,k) = cz_h_1(k)*kappaY_e(j)*dy*codt...
                                *murym1(i,j,k)/(kappaX_h(i)*dx);
          end
       end
    end
end


if(npmlz_2 > 0)
    b_hxy_zy_2 = zeros(nx-1,ny-1,npmlz_1);
    c_hx_zy_2 = zeros(nx-1,ny-1,npmlz_1);
    c_hy_zx_2 = zeros(nx-1,ny-1,npmlz_1);
    
    k = nz - npmlz_2;
    for k2 = 1:npmlz_2
       for j = 1:ny-1
          for i = 1:nx-1
             b_hxy_zy_2(i,j,k2) = bz_h_2(k2);
             c_hx_zy_2(i,j,k2) = cz_h_2(k2)*kappaX_e(i)*dx*codt...
                                 *murxm1(i,j,k)/(kappaY_h(j)*dy);
             c_hy_zx_2(i,j,k2) = cz_h_2(k2)*kappaY_e(j)*dy*codt...
                                 *murym1(i,j,k)/(kappaX_h(i)*dx);
          end
       end
       k = k+1;
    end
end

