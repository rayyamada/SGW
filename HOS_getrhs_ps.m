function R = HOS_getrhs_ps(S,zeta,kx,ky,g,h,M,t)
%
% - get RHS for HOS spectral code
%
%   eta_t = -eta_x PhiS_x - eta_y PhiS_y  + ( 1 + eta_x^2 + eta_y^2 ) Phi_z 
%  PhiS_t = -g eta - 0.5(PhiS_x^2 + PhiS_y^2) + 0.5(1 + eta_x^2 + eta_y^2) Phi_z^2  
%
% INPUT
%    S(:,:,1) = eta   (ps)
%    S(:,:,2) = PhiS  (ps)
%   eta - 2D field for surface elevation
%  PhiS - 2D field for surface potential evaluated at eta 
%  zeta - bottom undulation  (ps)
%
% kx,ky - wavenumber matricies
%     g - gravitational constant
%     h - mean domain height
%     M - perturbation order
%     t - time (not used)
%
% OUTPUT
%   R(:,:,1) = Reta
%   R(:,:,2) = RPhiS
%   Reta - rhs for  eta equation  (ps)
%  RPhiS - rhs for PhiS equation  (ps)


R = zeros(size(S));

 eta_ps = S(:,:,1);
PhiS_ps = S(:,:,2);

 eta = ps2k(eta_ps);
PhiS = ps2k(PhiS_ps);

% compute spectral derivatives
 eta_x = 1i*kx.*eta;
PhiS_x = 1i*kx.*PhiS; 

 eta_y = 1i*ky.*eta;
PhiS_y = 1i*ky.*PhiS;

%----------------------------------------------------------
% get Phi_z using method in Liu and Yue 1998
 Phi_z_ps = HOS_getPhiz(kx,ky,eta_ps,PhiS_ps,zeta,h,M);
%----------------------------------------------------------


 % convert spectral quantities to physical space (ps)
 eta_x_ps = k2ps(eta_x);
 eta_y_ps = k2ps(eta_y);
PhiS_x_ps = k2ps(PhiS_x);
PhiS_y_ps = k2ps(PhiS_y);
 
termA_ps = 1 + eta_x_ps.^2 + eta_y_ps.^2;  
termB_ps = termA_ps.*Phi_z_ps;
termC_ps = termA_ps.*(Phi_z_ps.^2);
termD_ps = PhiS_x_ps.^2 + PhiS_y_ps.^2;

% -- full nonlinear RHS -------------------------------
 Reta =  -eta_x_ps.*PhiS_x_ps - eta_y_ps.*PhiS_y_ps + termB_ps ;
RPhiS = -g*eta_ps  - 0.5*termD_ps + 0.5*termC_ps ;

% -- linear RHS ---------------------------------------
%  Reta =   Phi_z_ps ;
%  RPhiS = -g*eta_ps ;   
 
R(:,:,1) = Reta;
R(:,:,2) = RPhiS;

return
   
 
