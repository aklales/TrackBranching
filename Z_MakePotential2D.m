%% Make Morse Oscillator potential and initial wave function
% DESCRIPTIVE TEXT
clear all;


hbar = 1;



%%  coordinates
dx = 0.05;              % Real length between points on the grid
dy = 0.05;
nx = 4096;              % Number of grid points
ny = 1024;



xx = dx:dx:(dx*nx);
yy = dy:dy:(dy*ny);
Lx = max(xx);
Ly = max(yy);
M = length(xx);
N = length(yy);


%%  potential
sigma_rand = 1.5;
sigma_QPC = 3;
sigma_gap = 1;  % = 1, b = 128 corresponds to energy of 8
b = 128;

                                         %number of gaussians to do

[Vrand,V_QPC,GaussianCenters, Fx_rand, Fy_rand, Fxx_rand, Fyy_rand, Fxy_rand , Fx_QPC, Fy_QPC, Fxx_QPC, Fyy_QPC, Fxy_QPC]  = Random_Potential_with_QPC(xx,yy,5000,sigma_rand,b,sigma_QPC,sigma_gap);
Vimag = Z_AbsorbingBoundary(M,N,0.05);



V = Vrand + V_QPC + Vimag;


%%  Output information needed for Classical Calculation
save('Potential_4096_1024_i.mat','Vrand','V_QPC','Vimag','xx','yy', 'Fx_rand', 'Fy_rand', 'Fxx_rand', 'Fyy_rand', 'Fxy_rand' , 'Fx_QPC', 'Fy_QPC', 'Fxx_QPC', 'Fyy_QPC', 'Fxy_QPC');


