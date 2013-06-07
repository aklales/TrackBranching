function[Vrand,V_QPC,GaussianCenters, Fx_rand, Fy_rand, Fxx_rand, Fyy_rand, Fxy_rand , Fx_QPC, Fy_QPC, Fxx_QPC, Fyy_QPC, Fxy_QPC] = Random_Potential_with_QPC(xx,yy,nGauss,sigma_rand,b,sigma_QPC,sigma_gap)

%%  This function creates a random potential with nGauss gaussians as its base.  
%  For no QPC, set b = 0
%  xx and yy should be the arrays of positions along the width and the length of the sample
%%



m = length(xx);
n = length(yy);

Vrand = zeros(m,n);
Fx_rand = zeros(m,n);
Fy_rand = zeros(m,n);
Fxx_rand = zeros(m,n);
Fxy_rand = zeros(m,n);
Fyy_rand = zeros(m,n);
V_QPC = zeros(m,n);
Fx_QPC = zeros(m,n);
Fy_QPC = zeros(m,n);
Fxx_QPC = zeros(m,n);
Fxy_QPC = zeros(m,n);
Fyy_QPC = zeros(m,n);

%Force = zeros(5,m,n);



Lx = max(xx);
Ly = max(yy);


GaussianCenters = zeros(3,nGauss);

for ii = 1:nGauss
    GaussianCenters(1,ii) = Lx*rand(); 
    GaussianCenters(2,ii) = Ly*rand();
    GaussianCenters(3,ii) = 2*rand() - 1;
end

% For each point on the grid...
for ii = 1:m
    if mod(ii,50)==0
        ii
    end
    for jj = 1:n
        
        % We get effects from the Gaussians...
        [Fx,Fy,Fxx,Fyy,Fxy,pot] = Find__Force(xx(ii),yy(jj),GaussianCenters,sigma_rand);
        
        V_QPC(ii,jj) = b*(exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)) - exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2)));
        Fx_QPC(ii,jj) = b*(  (1/sigma_QPC^2)  *  (xx(ii)-2*sigma_QPC)  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2))  -  (1/sigma_QPC^2)  *  (xx(ii)-2*sigma_QPC)  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2)));
        Fy_QPC(ii,jj) = b*(- (1/sigma_gap^2)  *  (yy(jj)-Ly/2)  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2)));
        
        Fxx_QPC(ii,jj) = b*(  (1/sigma_QPC^2)^2  *  (xx(ii)-2*sigma_QPC)^2  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2))  -  (1/sigma_QPC^2)^2  *  (xx(ii)-2*sigma_QPC)^2  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2))  - (1/sigma_QPC^2)  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2))  +  (1/sigma_QPC^2)  *   exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2)));
        Fyy_QPC(ii,jj) = b*(- (1/sigma_gap^2)^2  *  (yy(jj)-Ly/2)^2  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2))  +  (1/sigma_gap^2)    *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2)));
        Fxy_QPC(ii,jj) = b*(- (1/sigma_gap^2)  *  (yy(jj)-Ly/2)  *  (1/sigma_QPC^2)  *  (xx(ii)-2*sigma_QPC)  *  exp(-(xx(ii)-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy(jj)-Ly/2)^2 /(2*sigma_gap^2)));
        
        
        
        Vrand(ii,jj) = pot;
        Fx_rand(ii,jj) = Fx;
        Fy_rand(ii,jj) = Fy;
        Fxx_rand(ii,jj) = Fxx;
        Fxy_rand(ii,jj) = Fxy;
        Fyy_rand(ii,jj) = Fyy;
        %%Force(:,ii,jj) = [Fx,Fy,Fxx,Fyy,Fxy];
        %%Need to add force from qpc to this value for classical
        %%simulations
    end
end

