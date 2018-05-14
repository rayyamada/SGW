% HOS_spectral_ps - main code to integrate equations with HOS method forward in time
%		  - must initialize variables before calling this script
%
% Pseudospectral code, High-Order-Spectral (HOS) method (Liu and Yue 1998), 
%   integrates nonlinear free-surface B.C.s (Zakharov 1968) in physical space, 
%   RK4 in time. Square with periodic B.C.s
% 
%   eta_t = -eta_x PhiS_x - eta_y PhiS_y  + ( 1 + eta_x^2 + eta_y^2 ) Phi_z 
%
%   PhiS_t = -g eta - 0.5(PhiS_x^2 + PhiS_y^2) + 0.5(1 + eta_x^2 + eta_y^2) Phi_z^2  
% 
% - subscripts (_x,_y,_z,_t) denote partial derivatives
% -  eta is the surface elevation
% - PhiS is the surface potential 
% - Phi_z is evaluated at the surface elevation eta [Phi_z = Phi_z(x,y,eta,t)]
%
% 

%%
Sout  = zeros(N,N,2,nframes);
% Phiz = zeros(N,N,nframes);      % (optional) store Phiz at z=eta for each timestep

S = zeros(N,N,2);
S(:,:,1) = eta_0;
S(:,:,2) = PhiS_0;

%
frame = 1;
% save initial condition to Sout
Sout(:,:,1,frame) =  eta_0;
Sout(:,:,2,frame) = PhiS_0;


%% Time step forward
disp('begun time integration ...')
for i = 1:nt
    disp(i)
    
       % (optional) store Phiz at z=eta
%        Phiz(:,:,frame) = HOS_getPhiz(kx,ky, S(:,:,1),S(:,:,2), zeta,h,M);
                

%-------- TIMESTEP AND DIFFUSE -----------

        %intermediate steps
        % Y1 = S;
        Rk  = HOS_getrhs_ps(S,zeta,kx,ky,g,h,M,t(i));
        
        Y2 = S + 0.5*dt*Rk;
        Rk2 = HOS_getrhs_ps(Y2,zeta,kx,ky,g,h,M,t(i)+0.5*dt);
        
        Y3 = S + 0.5*dt*Rk2;
        Rk3 = HOS_getrhs_ps(Y3,zeta,kx,ky,g,h,M,t(i)+0.5*dt);
        
        Y4 = S + dt*Rk3;
        Rk4 = HOS_getrhs_ps(Y4,zeta,kx,ky,g,h,M,t(i)+dt);
        
        % time step
        S = S + dt/6*(Rk + 2*Rk2 + 2*Rk3 + Rk4);
    
    
    % diffuse
    S(:,:,1) = k2ps( ps2k(S(:,:,1))./M2 );
    S(:,:,2) = k2ps( ps2k(S(:,:,2))./M2 );
    

   % ---- SAVE STEP ----
   if mod(i,savestep)==0
%        disp(i)
       frame = frame + 1;
       
       % save to Sout
       Sout(:,:,1,frame) =  S(:,:,1);
       Sout(:,:,2,frame) =  S(:,:,2);
       
   end

   
end % for i

%% 


