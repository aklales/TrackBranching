%% 2D FFT split operator code


%% parameters for random potential
clear all;
DoYouWantAMovie = 0;  %Set to 1 to create and save an animation. Set to 0 save time.
hbar = 1;


dt = .04;
time = 36;

PEHO = 8.0;  %potential energy from harmonic oscillator


load('Potential_4096_1024_i.mat','xx','yy','Vimag','V_QPC','Vrand')  % loads xx, yy and V
disp('Potential loaded...')

addpath /MatlabFunctions
Vimag = AbsorbingBoundary(length(xx),length(yy),0.08);
rmpath /MatlabFunctions

V = V_QPC -  0.0033*1i*Vimag + Vrand;

%V = -10*1i*Vnothing;
%V = 0.5 * Vharm;
%V(:,:) = 0.0;

M = length(xx);
N = length(yy);
dx = xx(2)-xx(1);
dy = yy(2)-yy(1);
% Lx = xx(end);
% Ly = yy(end)
Lx = xx(end)-xx(1);
Ly = yy(end)-yy(1);
px = -pi/dx:2*pi/M/dx:(pi/dx-2*pi/M);
py = -pi/dy:2*pi/N/dy:(pi/dy-2*pi/N);
[py,px] = meshgrid(py,px);
psquare = px.^2+py.^2;
K = exp(-1i*psquare*dt/(2*hbar));
Kshift = fftshift(K);








%%  Initial wavefunction conditions
B = 20;         %   Initial wavefunction amplitude
a = 1.5;          %   Initial wavefunction width
x0 = 6.0;       %   Initial center position of inital wavefunction
y0 = Ly/2;
p0 = 4.0;       %   Initial wavefunction momentum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eigen_E = (p0^2)/2 + PEHO;



fprintf('Dimensions: %f by %f (%d by %d)\n',Lx,Ly,M,N);
%speed of wave out of qpc
speed = sqrt(p0^2 + 2*PEHO);
fprintf('Wave speed out of QPC is %f\n',speed);
fprintf('Total time is %d\n',time);


vcontours = 0:1*B:20*B;
%vcontours = 1:20:400;



if DoYouWantAMovie == 0;
    disp('no movie')
end

%%  Define the initial wave function
disp('Creating initial Gaussian wave function...')
Psi = zeros(M,N);
for ii = 1:M
    %str = fprintf('%d of %d\n', ii, M);
    for jj = 1:N

        rrsquare = ((yy(jj)-y0)^2 + (xx(ii)-x0)^2 );
        Psi(ii,jj) = B*(exp( 1i*p0*xx(ii)-rrsquare/(2*a^2)));

    end
end    

Eigen = Psi;
%Eigen = abs(Psi).^2;
V = exp(-1i*V*dt/hbar);

if DoYouWantAMovie == 1
    aviobj = avifile('harmonic.avi','compression','None');
    fig = figure;
end


disp('Initializing time loop...')

for kk = 1:(time/dt)
    fprintf('Timestep %d of %d\n',kk,time/dt);

    Psi = ifft2(Kshift.*fft2(V.*Psi));
    Eigen = Eigen + (Psi * exp(1i*Eigen_E*kk*dt));
    %Eigen = Eigen + abs(Psi).^2;
    abs2 = abs(Psi).^2;
    %prob = sum(sum(abs2))
    if DoYouWantAMovie ==1
        contour(abs(Psi).^2,vcontours);
        axis([460 560 460 560]);
        F = getframe(fig);
        aviobj = addframe(aviobj,F);
    end  

end

save('QuantumOutput2D.mat','Eigen','Psi')

disp('Done!')
if DoYouWantAMovie == 1
    aviobj = close(aviobj);
end