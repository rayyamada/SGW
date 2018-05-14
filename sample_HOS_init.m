%% sample initialization script to run HOS

% ------- grid setup -------------
a = -125;  % domain is [a,b]x[a,b]
b =  125;
L =  b-a;  % domain length

N = 512;  % resolution

dx = L/N;

% define wavenumbers
k = 2*pi/L*[0:N/2, -N/2+1:-1];   
[kx, ky]=  meshgrid(k, k);  % spectral grid kx,ky
K = sqrt(kx.^2+ky.^2);      % wavevector magnitude

x_ = a + dx*(0:N-1);
[x,y] = meshgrid(x_, x_);   % spatial grid x,y

% -------- time --------------------   
dt = 0.25;       % time step
T  = 300;        % final time

 t = 0:dt:T;     % time for integration
nt = numel(t);

savestep = 1;    % save only when t = savestep*dt  
time = 0:savestep*dt:T;   % time variable corresponding to saved output
nframes = numel(time);

% -------- other parameters --------------------
% hyperdiffusion:  nu Laplacian^r u
nu = 1e-13;
r=4;          
  
M_disip = nu*dt *K.^(2*r) ;
M2 =  ones(N) + M_disip   ;  % used for hyperdiffusion computation

h = 1;       % mean domain height  (non-dimensionalized)
g = 1;       % gravity (non-dimensionalized)


% --- HOS --------------------------
M = 4;      % perturbation order


%% 1. initial conditions
% in this example we are creating a concave mirror topography setup
theta = atan(y./x);

lambda = L/40;          % 40 waves will fit in the domain
k1 = [2*pi/lambda, 0];  % initial plane wave

k1x = k1(1)*ones(N,N);
k1y = k1(2)*ones(N,N);

kbx = k1(1)*(1+cos(theta));  % topographic wavenumber kb
kby = k1(1)*sin(theta);

km = norm(k1,2);
w = sqrt(g*km*tanh(km*h));   % omega - frequency

P = 2*pi/w;       % period
c   = w/km;       % phase speed
cg  = g/(2*w)*( tanh(km*h) + km*h*(sech(km*h))^2 );    % group velocity dw/dkm


% initial wave amplitude  
A10 = 1e-2;

% topography height
d = 0.15;  

% ---topography
xb1 =  0*lambda;   % define boundaries of topographic patch [xb1,xb2]x[yb1,yb2] 
xb2 =  5*lambda;
yb1 = -5*lambda;  
yb2 =  5*lambda;


temp = sin(k1(1)*(sqrt(x.^2 + y.^2) + x));
temp( (x<xb1)|(x>xb2)|(y<yb1)|(y>yb2) ) = 0; 

zeta = d*temp;    % zeta is the topography with amplitude d   

% HOS initial conditions
 eta_0 =       -A10*sin(k1(1)*x + k1(2)*y) ;
PhiS_0 =    g/w*A10*cos(k1(1)*x + k1(2)*y) ;

%% 

% to run models the models do as follows:
%
% HOS_spectral_ps
%    
%  eta = squeeze(Sout(:,:,1,1:end-1));
% PhiS = squeeze(Sout(:,:,2,1:end-1));
%  

