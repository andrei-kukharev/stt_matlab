%{ 
The study of the STT-effect using the Landau–Lifshitz–Gilbert-Slonczewski-Berger equation.
%}

function tst

% matherial parameters
a = 0.02;  
b = 2.856; 
ha = 0.42;  
na = [0 1 0]; 
Nx = 0.32; Nz = 0.36; Ny=Nx;
Y0 = [11*pi/24; pi/24];  
J0 = 3.165995e-03; 
w1 = 0.1;
Nr = [Nx Ny Nz];
kf = 400;  
wr_in_file = 1;
Tm = 2000;  
nx = na(1); ny = na(2); nz = na(3); 
sx = 1; sy = 0; sz = 0;
kt = 3.2321e-12; 

hex = 0; hey = 0; hez = 0; % external field
f0 = 0.019; 
Nzp0 = 0; NzpA = 0.00; fNz =  0.8*f0;
Nzt = @(t) (Nz + Nzp0 + NzpA*sin(2*pi*fNz*t));

J0 = 0.12;
J = @(t) J0*((sign(sin(2*pi*f*t))+1)/2); % input current

% The magnetization dynamics equation
an = @(q,p) nx*sin(q)*cos(p) + ny*sin(q)*sin(p) + nz*cos(q);
g = @(q,p) 1/((3+ sx*sin(q)*cos(p)+ sy*sin(q)*sin(p) + sz*cos(q))*b - 4);

hx = @(q,p,t) hex + ha*nx*an(q,p) -(J(t)*g(q,p))*(sz*sin(q)*sin(p)-sy*cos(q)) - Nx*sin(q)*cos(p);
hy = @(q,p,t) hey + ha*ny*an(q,p) -(J(t)*g(q,p))*(sx*cos(q)-sz*sin(q)*cos(p)) - Ny*sin(q)*sin(p);
hz = @(q,p,t) hez + ha*nz*an(q,p) -(J(t)*g(q,p))*(sy*cos(p)-sx*sin(p))*sin(q) - Nzt(t)*cos(q);

hp = @(q,p,t) (-hx(q,p,t)*sin(p) + hy(q,p,t)*cos(p));
hq = @(q,p,t) hx(q,p,t)*cos(q)*cos(p) + hy(q,p,t)*cos(q)*sin(p) - hz(q,p,t)*sin(q);

osc = @(t,y) [ hp(y(1),y(2),t) + a*hq(y(1),y(2),t) ;
    (a*hp(y(1),y(2),t) - hq(y(1),y(2),t) )/sin(y(1)) ];    

options = odeset('MaxStep', 0.2);
[T,Y] = ode45(osc, [0 Tm], Y0, options); % ODE solution

% Converting from spherical coordinates to Cartesian coordinates
mx = sin(Y(:,1)).*cos(Y(:,2));
my = sin(Y(:,1)).*sin(Y(:,2));
mz = cos(Y(:,1));

Jt = J(T);
NztT = Nzt(T);
Jmax = max(abs(Jt));

% Fourier analysis:
n = 100000;  % the number of points
fx = fft(mx,n);
fy = fft(my,n);
fz = fft(mz,n);
fj = fft(Jt,n);
ppy = fy.* conj(fy) / n;
ppx = fx.* conj(fx) / n;
ppz = fz.* conj(fz) / n;
ppj = fj.* conj(fj) / n;

l=length(T);
dt=Tm/(l-1);
f = 1/dt*(0:n/kf)/n;  
px = ppx(1:n/kf+1);
py = ppy(1:n/kf+1);
pz = ppz(1:n/kf+1);
pj = ppj(1:n/kf+1);

% Visualization:
figure(1);

subplot(2,2,1);
plot(T,mx,'r-');
axis([0 Tm -1.01 1.01]);
xlabel('t');
ylabel('Mx');
grid on;

subplot(2,2,2);
plot(f,px)
grid on;
title('spectrum')
xlabel('f, a.i.')
ylabel('Ampl. mx');

subplot(2,2,3);
plot(T,NztT,'g-');
axis([0 Tm 0.00 3.00]);
title(strcat('Freq. f=',num2str(fNz)));
xlabel('t');
ylabel('Nz(t)');
grid on;

subplot(2,2,4);
plot3(mx,my,mz);
axis([-1.01 1.01 -1.01 1.01 -1.01 1.01]);
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
print -dpng -r100 'out\spectre1.png'; 

% Saving the result into a file
D_spx = [f/kt*1e-9; px.'];
fp = fopen('out\spectre_x_GHz.txt','wt');  
fprintf(fp,'%f %f\n',D_spx);
fclose(fp);

end

