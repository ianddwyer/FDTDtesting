%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialize the fields (Default is zero)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%
%  E-Fields
%%%%%%%%%%%%%%%%

ex = zeros(nx-1,ny,nz);
ey = zeros(nx,ny-1,nz);
ez = zeros(nx,ny,nz-1);

%%%%%%%%%%%%%%%%
%  H-Fields
%%%%%%%%%%%%%%%%

hx = zeros(nx,ny-1,nz-1);
hy = zeros(nx-1,ny,nz-1);
hz = zeros(nx-1,ny-1,nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PML Auxiliary Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%  PML E-Field-updates
%%%%%%%%%%%%%%%%

if(npmlx_1 > 0)
    q_ey_xz_1 = zeros(npmlx_1+1,ny-1,nz-1);
    q_ez_xy_1 = zeros(npmlx_1+1,ny-1,nz-1);
end

if(npmlx_2 > 0)
    q_ey_xz_2 = zeros(npmlx_2,ny-1,nz-1);
    q_ez_xy_2 = zeros(npmlx_2,ny-1,nz-1);
end

if(npmly_1 > 0)
    q_ez_yx_1 = zeros(nx-1,npmly_1+1,nz-1);
    q_ex_yz_1 = zeros(nx-1,npmly_1+1,nz-1);
end

if(npmly_2 > 0)
    q_ez_yx_2 = zeros(nx-1,npmly_2,nz-1);
    q_ex_yz_2 = zeros(nx-1,npmly_2,nz-1);
end


if(npmlz_1 > 0)
    q_ex_zy_1 = zeros(nx-1,ny-1,npmlz_1+1);
    q_ey_zx_1 = zeros(nx-1,ny-1,npmlz_1+1);
end

if(npmlz_2 > 0)
    q_ex_zy_2 = zeros(nx-1,ny-1,npmlz_2);
    q_ey_zx_2 = zeros(nx-1,ny-1,npmlz_2);
end


%%%%%%%%%%%%%%%%
%  PML H-Field-updates
%%%%%%%%%%%%%%%%

if(npmlx_1 > 0)
    q_hy_xz_1 = zeros(npmlx_1,ny-1,nz-1);
    q_hz_xy_1 = zeros(npmlx_1,ny-1,nz-1);
end

if(npmlx_2 > 0)
    q_hy_xz_2 = zeros(npmlx_2,ny-1,nz-1);
    q_hz_xy_2 = zeros(npmlx_2,ny-1,nz-1);
end


if(npmly_1 > 0)
    q_hz_yx_1 = zeros(nx-1,npmly_1,nz-1);
    q_hx_yz_1 = zeros(nx-1,npmly_1,nz-1);
end

if(npmly_2 > 0)
    q_hz_yx_2 = zeros(nx-1,npmly_2,nz-1);
    q_hx_yz_2 = zeros(nx-1,npmly_2,nz-1);
end

if(npmlz_1 > 0)
    q_hx_zy_1 = zeros(nx-1,ny-1,npmlz_1);
    q_hy_zx_1 = zeros(nx-1,ny-1,npmlz_1);
end

if(npmlz_2 > 0)
    q_hx_zy_2 = zeros(nx-1,ny-1,npmlz_2);
    q_hy_zx_2 = zeros(nx-1,ny-1,npmlz_2);
end
