function [Phi_z, PhiSm] = HOS_getPhiz(kx,ky,eta,PhiS,zeta,h,M)

% INPUTs
%   kx,ky = wavenumber domain
%
%   eta =  eta(x,y)  at time t   - surface elevation (ps)
%  PhiS = PhiS(x,y)  at time t   - surface potential (ps)
%  zeta = zeta(x,y)              - bottom undulation (ps)
%     h = mean height of domain
%     M = perturbation order
%
% OUTPUT
%   Phi_z = Phi_z(x,y) at time t (ps)
%   PhiSm = PhiSm(x,y,M) at time t,  stores PhiS^(m) for orders m=1,...,M
%

N = size(kx,1);  % domain is square, NxN

% Phi_z = zeros(N,N);   % (x,y)

K = sqrt(kx.^2 + ky.^2);

% quantities used for calculation
alpha_ = zeros(N,N,M);   % alpha_n1,n2^(m) quantity needed to compute (5.5) 
 beta_ = zeros(N,N,M);   %  beta_n1,n2^(m) quantity needed to compute (5.6)   
 
    
    dzl_Phi_m = zeros(N,N,M);  
        PhiSm = zeros(N,N,M);    % store PhiS^(m) for each order m

%--------------------------------------------------------------                                
%--------------------------------------------------------------                                
% --- GET alpha_n^(m)  and  beta_n^(m) ---    

% m = 1
       alpha_(:,:,1) = 1/N^2 * fft2( PhiS );
%       beta_(:,:,1) = 0; 
       
% m = 2,3,...,M
for m = 2:M
    termC  = zeros(N,N);  % (x,y)
    for l = 1:m-1
        
        dzlA = zeros(N,N);  % (x,y)
        dzlB = zeros(N,N);
        
        if mod(l,2) == 1      % ----- l is odd
            dzlA = real( N^2*ifft2( alpha_(:,:,m-l).*(K.^l).*tanh(K*h) ));
            dzlB = real( N^2*ifft2(  beta_(:,:,m-l).*(K.^(l-1)).*(1./cosh(K*h)) )); 
        elseif mod(l,2) == 0  % ----- l is even
            dzlA = real( N^2*ifft2( alpha_(:,:,m-l).*(K.^l) ));
          % dzlB = zeros(N,N); 
        end
          
        % compute sum over l
        termC = termC + eta.^(l)/factorial(l).*(dzlA+dzlB);
      
    end % l
        
     alpha_(:,:,m) = -1/N^2 * fft2( termC );
     
   %% --------------------------
  
     termC2 = zeros(N,N);   
     Kn = K;
     Kn(1,1) = 1;  % to prevent 1/0 = NaN
                   % this term will be zeroed out anyway
     for l_ = 1:m-1
        l = l_ - 1;   % l from 0 to m-2
                      % l is the order of the lth derivative

        dzlA_x = zeros(N,N);   % (x,y)
        dzlB_x = zeros(N,N);   
        dzlA_y = zeros(N,N);   
        dzlB_y = zeros(N,N);   
             
        % dzlA_x,y; dzlB_x,y
        if mod(l,2) == 1     % ----- l is odd
          % dzlA_x = zeros(N,N);
          % dzlA_y = zeros(N,N);
            
            dzlB_x =  real(N^2*ifft2( ...
                (1i*kx).* beta_(:,:,m-l_).*(K.^(l-1)) ));
            dzlB_y =  real(N^2*ifft2( ...
                (1i*ky).* beta_(:,:,m-l_).*(K.^(l-1)) ));
        elseif mod(l,2) == 0 % ----- l is even
            dzlA_x =  real(N^2*ifft2( ...
                (1i*kx).*alpha_(:,:,m-l_).*(K.^l).*(1./cosh(K*h)) ));
            dzlA_y =  real(N^2*ifft2( ...
                (1i*ky).*alpha_(:,:,m-l_).*(K.^l).*(1./cosh(K*h)) ));
            
            dzlB_x = -real(N^2*ifft2( ...
                (1i*kx).* beta_(:,:,m-l_).*(Kn.^(l-1)).*tanh(K*h) ));
            dzlB_y = -real(N^2*ifft2( ...
                (1i*ky).* beta_(:,:,m-l_).*(Kn.^(l-1)).*tanh(K*h) ));
           % need to use Kn here because if l=0, will get 1/0 when 
           % (n1,n2) = (0,0).  But this term will be zeroed out anyway
           % from tanh(0) and (i1*0) terms.
           %
           % note: ifft2 sums over (n1,n2)=(0,0), but because of derivative
           % these terms are all zero from the i1*kx, 1i*ky terms in front
        end
              
      dzl_Phi_x = dzlA_x + dzlB_x ;
      dzl_Phi_y = dzlA_y + dzlB_y ;
      
      termA2 = zeta.^l_ /factorial(l_) .* dzl_Phi_x;  % product in ps
      termB2 = zeta.^l_ /factorial(l_) .* dzl_Phi_y;   
                                                
      termA2_x = k2ps( 1i*kx.*ps2k( termA2 ) ); % x derivative of termA
      termB2_y = k2ps( 1i*ky.*ps2k( termB2 ) ); % y derivative of termB
      
      termC2 = termC2 + (termA2_x + termB2_y);   % summation over l_    
      
     end % l_
     
     beta_(:,:,m) = 1/N^2 * fft2( termC2 );
   %--------------------------
  
end   % m


%--------------------------------------------------------------
%--------------------------------------------------------------
% --- EQN (5.7) --- get dzl_Phi_m
% using different summation method which separates orders clearly
for m = 1:M
    
    for l = 0:m-1
        lp1 = l+1;
        
        dzlA = zeros(N,N);  % use for lth derivative of Phi^m
        dzlB = zeros(N,N);
        
        dzlp1A = zeros(N,N); % use for (l+1)th derivative of Phi^m
        dzlp1B = zeros(N,N);
        
        if mod(l,2) == 1      % ----- l is odd
            dzlA = real( N^2*ifft2( alpha_(:,:,m-l).*(K.^l).*tanh(K*h) ));
            dzlB = real( N^2*ifft2(  beta_(:,:,m-l).*(K.^(l-1)).*(1./cosh(K*h)) )); 
        elseif mod(l,2) == 0  % ----- l is even
            dzlA = real( N^2*ifft2( alpha_(:,:,m-l).*(K.^l) ));
          % dzlB = zeros(N,N); 
        end
        
        if mod(lp1,2) == 1      % ----- l is odd
            dzlp1A = real( N^2*ifft2( alpha_(:,:,m-l).*(K.^lp1).*tanh(K*h) ));
            dzlp1B = real( N^2*ifft2(  beta_(:,:,m-l).*(K.^(lp1-1)).*(1./cosh(K*h)) )); 
        elseif mod(lp1,2) == 0  % ----- l is even
            dzlp1A = real( N^2*ifft2( alpha_(:,:,m-l).*(K.^lp1) ));
          % dzlp1B = zeros(N,N); 
        end
        
        % multiply by factor of eta^l_/factorial(l_) in (5.7)
        temp = eta.^(l)/factorial(l);
        
        
            PhiSm(:,:,m) =     PhiSm(:,:,m) + temp.*(dzlA + dzlB);
        dzl_Phi_m(:,:,m) = dzl_Phi_m(:,:,m) + temp.*(dzlp1A + dzlp1B);
        
    end % l
        
end % m

Phi_z  = sum(dzl_Phi_m,3);   


% END




