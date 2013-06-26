%% Make Morse Oscillator potential and initial wave function
%  This code makes a potential grid to feed into the Quantum Simulation and
%  makes a list of gaussians and their centers to feed into the Classical
%  Simulation

clear all;
WriteProgress = 1;

%%  coordinates
dx = 0.05;              % Real length between points on the grid
dy = 0.05;
nx = 4096;              % Number of grid points
ny = 2048;

%%  Potential Characteristics
sigma_rand = 1.5;
sigma_QPC = 3;
sigma_gap = 1;  % = 1, b = 128 corresponds to energy of 8
b = 128;
nGauss = 10000;



%%Define major output objects
xx = dx:dx:(dx*nx);
yy = dy:dy:(dy*ny);
GaussianCenters = zeros(3,nGauss);
Potential = zeros(nx,ny);




Lx = max(xx);
Ly = max(yy);
M = length(xx);
N = length(yy);

if WriteProgress == 1
    fprintf('Calculating %d Gaussians\n',nGauss);
end


%% Define Potential
for ii = 1:nGauss
    
    GaussianCenters(1,ii) = Lx*rand(); 
    GaussianCenters(2,ii) = Ly*rand();
    GaussianCenters(3,ii) = 2*rand() - 1;
end

%%  Calculate potential grid
if WriteProgress == 1
    fprintf('Calculating Potential Grid, %d by %d\n',nx,ny);
end

addpath ./MatlabFunctions
for ii = 1:nx
    if (WriteProgress == 1 && mod(ii,10)==0)
        fprintf('%d\n',ii);
    end
    
    
    for jj = 1:ny
        Potential(ii,jj) = LocalForce(xx(ii),yy(jj),GaussianCenters,sigma_rand,b,sigma_QPC,sigma_gap,Ly,0);
    end
end
rmpath ./MatlabFunctions

%Vimag = Z_AbsorbingBoundary(M,N,0.05);
%V = Vrand + V_QPC + Vimag;


%%  Output information needed for Classical Calculation
%save('Potential_4096_1024_i.mat','GaussianCenters','Vrand','V_QPC','Vimag','xx','yy', 'Fx_rand', 'Fy_rand', 'Fxx_rand', 'Fyy_rand', 'Fxy_rand' , 'Fx_QPC', 'Fy_QPC', 'Fxx_QPC', 'Fyy_QPC', 'Fxy_QPC','sigma_rand','b','sigma_gap','sigma_QPC');

save('Potential_4096_2048_TEST.mat','GaussianCenters','Potential','xx','yy','sigma_rand','b','sigma_gap','sigma_QPC');
