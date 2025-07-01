%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize PML update parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E-Update coefficients (Ampere's law):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X-directed PML:
kappaX_e = ones(nx,1);
kappaX_h = ones(nx-1,1);

% X-Minus boundary:
%%%%%%%%%%%%%%%%%%%%
del = dx;
iend = 1;
npml = npmlXMinus;  % PML Depth in cells
npmlx_1 = npml;
m = mScaleXMinus;     % mscaling of sigma, kappa
sigMxScal = pSigmaOptXMinus;  % scaling of sigma-optimal
kapmax = kappaMaxXMinus;     % kappa-max
alpmax = aMaxXMinus;     % a-max
ma = aScalXMinus;     % a scale
epsrEff = epsrEffXMinus;     % epsilon effective

% estimate sigma-optimal, and then calculate sigma-Max
sigOpt = (m+1.0)/(150*pi*sqrt(epsrEff)*del);
sigMax = sigMxScal*sigOpt;

[kappae, bx_e_1, cx_e_1] = ...
  PmlScaleEUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaX_e(1:npml+1) = kappae;

[kappah, bx_h_1, cx_h_1] = ...
  PmlScaleHUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaX_h(1:npml) = kappah;
clear kappae;
clear kappah;

% X-Plus boundary:
%%%%%%%%%%%%%%%%%%%%
iend = 2;
npml = npmlXPlus;  % PML Depth in cells
npmlx_2 = npml;
m = mScaleXPlus;     % mscaling of sigma, kappa
sigMxScal = pSigmaOptXPlus;  % scaling of sigma-optimal
kapmax = kappaMaxXPlus;     % kappa-max
alpmax = aMaxXPlus;     % a-max
ma = aScalXPlus;     % a scale
epsrEff = epsrEffXPlus;     % epsilon effective

% estimate sigma-optimal, and then calculate sigma-Max
sigOpt = (m+1.0)/(150*pi*sqrt(epsrEff)*del);
sigMax = sigMxScal*sigOpt;

[kappae, bx_e_2, cx_e_2] = ...
  PmlScaleEUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaX_e(nx-npml:nx) = kappae;

[kappah, bx_h_2, cx_h_2] = ...
  PmlScaleHUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaX_h(nx-npml:nx-1) = kappah;
clear kappae;
clear kappah;



% Y-directed PML:
kappaY_e = ones(ny,1);
kappaY_h = ones(ny-1,1);

% Y-Minus boundary:
%%%%%%%%%%%%%%%%%%%%
del = dy;
iend = 1;
npml = npmlYMinus;  % PML Depth in cells
npmly_1 = npml;
m = mScaleYMinus;     % mscaling of sigma, kappa
sigMxScal = pSigmaOptYMinus;  % scaling of sigma-optimal
kapmax = kappaMaxYMinus;     % kappa-max
alpmax = aMaxYMinus;     % a-max
ma = aScalYMinus;     % a scale
epsrEff = epsrEffYMinus;     % epsilon effective

% estimate sigma-optimal, and then calculate sigma-Max
sigOpt = (m+1.0)/(150*pi*sqrt(epsrEff)*del);
sigMax = sigMxScal*sigOpt;

[kappae, by_e_1, cy_e_1] = ...
  PmlScaleEUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaY_e(1:npml+1) = kappae;

[kappah, by_h_1, cy_h_1] = ...
  PmlScaleHUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaY_h(1:npml) = kappah;


% Y-Plus boundary:
%%%%%%%%%%%%%%%%%%%%
iend = 2;
npml = npmlYPlus;  % PML Depth in cells
npmly_2 = npml;
m = mScaleYPlus;     % mscaling of sigma, kappa
sigMxScal = pSigmaOptYPlus;  % scaling of sigma-optimal
kapmax = kappaMaxYPlus;     % kappa-max
alpmax = aMaxYPlus;     % a-max
ma = aScalYPlus;     % a scale
epsrEff = epsrEffYPlus;     % epsilon effective

% estimate sigma-optimal, and then calculate sigma-Max
sigOpt = (m+1.0)/(150*pi*sqrt(epsrEff)*del);
sigMax = sigMxScal*sigOpt;

[kappae, by_e_2, cy_e_2] = ...
  PmlScaleEUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaY_e(ny-npml:ny) = kappae;

[kappah, by_h_2, cy_h_2] = ...
  PmlScaleHUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaY_h(ny-npml:ny-1) = kappah;
clear kappae;
clear kappah;



% Z-directed PML:
kappaZ_e = ones(nz,1);
kappaZ_h = ones(nz-1,1);

% Z-Minus boundary:
%%%%%%%%%%%%%%%%%%%%
del = dz;
iend = 1;
npml = npmlZMinus;  % PML Depth in cells
npmlz_1 = npml;
m = mScaleZMinus;     % mscaling of sigma, kappa
sigMxScal = pSigmaOptZMinus;  % scaling of sigma-optimal
kapmax = kappaMaxZMinus;     % kappa-max
alpmax = aMaxZMinus;     % a-max
ma = aScalZMinus;     % a scale
epsrEff = epsrEffZMinus;     % epsilon effective

% estimate sigma-optimal, and then calculate sigma-Max
sigOpt = (m+1.0)/(150*pi*sqrt(epsrEff)*del);
sigMax = sigMxScal*sigOpt;

[kappae, bz_e_1, cz_e_1] = ...
  PmlScaleEUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaZ_e(1:npml+1) = kappae;

[kappah, bz_h_1, cz_h_1] = ...
  PmlScaleHUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaZ_h(1:npml) = kappah;


% Z-Plus boundary:
%%%%%%%%%%%%%%%%%%%%
iend = 2;
npml = npmlZPlus;  % PML Depth in cells
npmlz_2 = npml;
m = mScaleZPlus;     % mscaling of sigma, kappa
sigMxScal = pSigmaOptZPlus;  % scaling of sigma-optimal
kapmax = kappaMaxZPlus;     % kappa-max
alpmax = aMaxZPlus;     % a-max
ma = aScalZPlus;     % a scale
epsrEff = epsrEffZPlus;     % epsilon effective

% estimate sigma-optimal, and then calculate sigma-Max
sigOpt = (m+1.0)/(150*pi*sqrt(epsrEff)*del);
sigMax = sigMxScal*sigOpt;

[kappa, bz_e_2, cz_e_2] = ...
  PmlScaleEUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaZ_e(nz-npml:nz) = kappa;

[kappah, bz_h_2, cz_h_2] = ...
  PmlScaleHUpdate(iend,del,dt,eps0,npml,m,ma,sigMax,kapmax,alpmax);

kappaZ_h(nz-npml:nz-1) = kappah;

