%%  Propagate the particles
%   This code classically propagates np particles over our potential
clear all;


DontReload = 1;  %if 0, reloads potential from file; if 1, generates potential on the fly from Gaussian Centers; if otherwise, uses potential already loaded into memory



if DontReload == 0
    %clear all;
    fprintf('Potential is loading...\n')
    load('Potential_4096_1024_i.mat')%loads Force and coords of entries in force matrix
    fprintf('Potential loaded.\n')
    
    nW = length(xx);
    nL = length(yy);
    xV = xx;
    yV = yy;
    Fx = Fx_rand + Fx_QPC;
    Fy = Fy_rand + Fy_QPC;
    V = Vrand + V_QPC;

    dx = abs(xV(2)-xV(1));  %spatial difference between adjacent points of potential
    dy = abs(yV(2)-yV(1));
    Lx = xV(end);           %real dimensions of rectangle 
    Ly = yV(end);


end

if DontReload ==1
    %clear all;
    fprintf('Gaussian Centers are loading...\n')
    load('Potential_4096_1024_i_GS.mat','xx','yy','GaussianCenters', 'sigma_rand', 'b','sigma_QPC','sigma_gap')%loads Force and coords of entries in force matrix
    fprintf('loaded!\n')
    
    nW = length(xx);
    nL = length(yy);
    xV = xx;
    yV = yy;


    dx = abs(xV(2)-xV(1));  %spatial difference between adjacent points of potential
    dy = abs(yV(2)-yV(1));
    Lx = xV(end);           %real dimensions of rectangle 
    Ly = yV(end);
    
end


D2 = 1;
m = 1;
DoStability = 1;
PEHO = 8;

np = 100;


vel = 4;                % initial velocity of particle
time = 25000;
dt = 0.001;

vf = sqrt(2*PEHO + vel^2);
buffx = ceil(vf*dt/dx) + 1;
buffy = ceil(vf*dt/dy) + 1;



if (DoStability == 1 && DontReload == 0)
    fprintf('Stability will be calculated.\n')
    Fxx = Fxx_rand + Fxx_QPC;
    Fyy = Fyy_rand + Fyy_QPC;
    Fxy = Fxy_rand + Fxy_QPC;
end

xx = zeros(np,time); 
yy = zeros(np,time);

CalculateVelocity = 1;
if CalculateVelocity == 1
    vvxx = zeros(np,time);
    vvyy = zeros(np,time);
end

%%  SET UP INTIAL CONDITIONS
fprintf('Initializing...\n')
    

vy = zeros(np,1);           %initial momenta of particles
vx =  vel*ones(np,1); 

%  (center location, standard deviation, number of particles, nothing)
xxi = normrnd(6,1.5,np,1);
yyi = normrnd(Ly/2,1.5,np,1);
if (DontReload ~= 1)
    for ii = 1:np
    	indx = round(xxi(ii)/dx);
        indy = round(yyi(ii)/dy);
        while indx < 1
            xxi(ii) = normrnd(6,1.5);
            indx = round(xxi(ii)/dx);
        end
        while indy < 1
            yyi(ii) = normrnd(Ly/2,1.5);
            indy = round(yyi(ii)/dy);
        end
        while (V(indx,indy) > PEHO)
            xxi(ii) = normrnd(6,1.5);
            yyi(ii) = normrnd(Ly/2,1.5);
            indx = round(xxi(ii)/dx);
            indy = round(yyi(ii)/dy);
            while indx < 1
                xxi(ii) = normrnd(6,1.5);
                indx = round(xxi(ii)/dx);
                bork = 1
            end
        end
        vy(ii,1) = sign(rand - 0.5) * sqrt((PEHO - V(indx,indy))/(2*m));
    
    end
else
    for ii = 1:np

        while xxi(ii) < dx
            xxi(ii) = normrnd(6,1.5);
        end
        while yyi(ii) < dy
            yyi(ii) = normrnd(Ly/2,1.5);
        end
        while (LocalForce(xxi(ii),yyi(ii),GaussianCenters,sigma_rand,b,sigma_QPC,sigma_gap,Ly,0) > PEHO)
            xxi(ii) = normrnd(6,1.5);
            yyi(ii) = normrnd(Ly/2,1.5);
            while xxi(ii) < dx
                xxi(ii) = normrnd(6,1.5);
                bork = 1
            end
        end
        vy(ii,1) = sign(rand - 0.5) * sqrt((PEHO - LocalForce(xxi(ii),yyi(ii),GaussianCenters,sigma_rand,b,sigma_QPC,sigma_gap,Ly,0))/(2*m));
    
    end
end


if CalculateVelocity == 1
    vvxx(:,1) = vx;
    vvyy(:,1) = vy;
end


xx(:,1) = xxi;
yy(:,1) = yyi;

if DoStability ==1
    M = zeros(np,time,4,4);
    M(:,:,1,1) = 1;
    M(:,:,2,2) = 1;
    M(:,:,3,3) = 1;
    M(:,:,4,4) = 1;
end

%% BEGIN PARTICLE-TIME LOOPS

if DontReload == 1
    %calculate potential on the fly
    fprintf('Propagation has begun. Calculating potential on the fly.\n')
    for ii = 1:np
        if mod(ii,1)==0
            fprintf('On particle %d of %d\n',ii,np)
        end
        for jj = 1:time

        
            xr = round(xx(ii,jj)/dx);
            yr = round(yy(ii,jj)/dy);

            if ((xr < buffx) || (xr > (nW-buffx)) || (yr < buffy) || (yr > (nL-buffy))) 
                
                break 
            end                
            
            forces = LocalForce(xx(ii,jj),yy(ii,jj),GaussianCenters,sigma_rand,b,sigma_QPC,sigma_gap,Ly,1);
            xforce = forces(1);
            yforce = forces(2);
        
            vx(ii) = vx(ii) + dt*xforce/m;
            vy(ii) = vy(ii) + dt*yforce/m; 
            xx(ii,jj+1) = xx(ii,jj) + vx(ii)*dt;
            yy(ii,jj+1) = yy(ii,jj) + vy(ii)*dt;  
        
            if DoStability == 1
                 
                curvatures = LocalForce(xx(ii,jj),yy(ii,jj),GaussianCenters,sigma_rand,b,sigma_QPC,sigma_gap,Ly,2);
                Fxx = curvatures(1);
                Fxy = curvatures(2);
                Fyy = curvatures(3);
                K = zeros(4);
                K(3,1) = 1;
                K(4,2) = 1;
                K(1,3) = -Fxx;
                K(1,4) = -Fxy;
                K(2,3) = -Fxy;
                K(2,4) = -Fyy;
            
                msub = squeeze(M(ii,jj,:,:));
                %K1 = K*msub;
        
            
            
            
%                 K2 = zeros(4);
%                 K(3,1) = 1;
%                 K(4,2) = 1;
%                 K(1,3) = -Fxx(xf,yf);
%                 K(1,4) = -Fxy(xf,yf);
%                 K(2,3) = -Fxy(xf,yf);
%                 K(2,4) = -Fyy(xf,yf);

                if mod(jj,750) == 0
                    current_determinant = det(msub)
                end
                
                M(ii,jj+1,:,:) = msub + K*msub*dt;
            end
        
        
            if CalculateVelocity == 1
                vvxx(ii,jj) = vx(ii);
                vvyy(ii,jj) = vy(ii);
            end
        
        

        


        end
    end  
else
    %Use potential from file
    fprintf('Propagation has begun. Using potential from file.\n')
    for ii = 1:np
        if mod(ii,1)==0
            fprintf('On particle %d of %d\n',ii,np)
        end
        for jj = 1:time

        
            xf = floor(xx(ii,jj)/dx);
            xc = ceil(xx(ii,jj)/dx);
            yf = floor(yy(ii,jj)/dy);
            yc = ceil(yy(ii,jj)/dy);
            xr = round(xx(ii,jj)/dx);
            yr = round(yy(ii,jj)/dy);

            if ((xf < buffx) || (xf > (nW-buffx)) || (yf < buffy) || (yf > (nL-buffy))) 
                
                break 
            end                
            xforce = ((xc-(xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fx(xf,yf) + (-1*xf + (xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fx(xc,yc) + (xc-(xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fx(xf,yc) + (-1*xf + (xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fx(xc,yf));
            yforce = ((xc-(xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fy(xf,yf) + (-1*xf + (xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fy(xc,yc) + (xc-(xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fy(xf,yc) + (-1*xf + (xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fy(xc,yf));
        
        
        
            vx(ii) = vx(ii) + dt*xforce/m;
            vy(ii) = vy(ii) + dt*yforce/m; 
            xx(ii,jj+1) = xx(ii,jj) + vx(ii)*dt;
            yy(ii,jj+1) = yy(ii,jj) + vy(ii)*dt;  
        
            if DoStability == 1
                %xxmid = (xx(ii,jj+1)+xx(ii,jj))/2;
                %yymid = (yy(ii,jj+1)+yy(ii,jj))/2;
            
                K = zeros(4);
                K(3,1) = 1;
                K(4,2) = 1;
                K(1,3) = -Fxx(xf,yf);
                K(1,4) = -Fxy(xf,yf);
                K(2,3) = -Fxy(xf,yf);
                K(2,4) = -Fyy(xf,yf);
            
                msub = squeeze(M(ii,jj,:,:));
                %K1 = K*msub;
        
            
            
            
%                 K2 = zeros(4);
%                 K(3,1) = 1;
%                 K(4,2) = 1;
%                 K(1,3) = -Fxx(xf,yf);
%                 K(1,4) = -Fxy(xf,yf);
%                 K(2,3) = -Fxy(xf,yf);
%                 K(2,4) = -Fyy(xf,yf);
                %det(msub)
                M(ii,jj+1,:,:) = msub + K*msub*dt;
            end
        
        
            if CalculateVelocity == 1
                vvxx(ii,jj) = vx(ii);
                vvyy(ii,jj) = vy(ii);
            end
        
        

        


        end
    end 
    
    
end

xx = transpose(xx);
yy = transpose(yy);


save('Results.mat','xx','yy');

density = betterbinz(1024,256,xx,yy,Lx,Ly);











