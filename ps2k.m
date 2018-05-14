function F = ps2k(F_ps)
% - transform field from physical space to spectral

F = fft2(F_ps);


% OLD VERSION WITH ZERO-PADDING FOR ANTI-ALIASING
% - zero padding for anti-aliasing, truncate matrix back to size N^2
%   from size (3N/2)^2
% 
% F = (2/3)^2 *fft2(F_ps);
% F = cat(1, ...
%            cat(2,     F(1:N/2,1:N/2),     F(1:N/2,N+1:3*N/2) ), ...
%            cat(2, F(N+1:3*N/2,1:N/2), F(N+1:3*N/2,N+1:3*N/2) ) ...
%        );
