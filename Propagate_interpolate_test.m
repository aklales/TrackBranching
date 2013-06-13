%%  Propagate the particles
%   This code classically propagates np particles over our potential



DontReload = 0;



if DontReload == 0
    clear all;
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


D2 = 1;
m = 1;
DoStability = 1;
PEHO = 8;

np = 10;


vel = 4;                % initial velocity of particle
time = 250000;
dt = 0.0001;
buffx = ceil(vel*dt/dx);
buffy = ceil(vel*dt/dy);



if DoStability == 1
    fprintf('Stability will be calculated.\n')
    Fxx = Fxx_rand + Fxx_QPC;
    Fyy = Fyy_rand + Fyy_QPC;
    Fxy = Fxy_rand + Fxy_QPC;
end

targ = 1;
xx = zeros(np,time/targ); 
yy = zeros(np,time/targ);

CalculateVelocity = 0;
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
fprintf('Propagation has begun.\n')
for ii = 1:np
    if mod(ii,1)==0
        fprintf('On particle %d of %d\n',ii,np)
    end
    for jj = 1:time

        
        xf = floor(xxi(ii)/dx);
        xc = ceil(xxi(ii)/dx);
        yf = floor(yyi(ii)/dy);
        yc = ceil(yyi(ii)/dy);

        if ((xf < buffx) || (xf > (nW-buffx)) || (yf < buffy) || (yf > (nL-buffy))) 
                
            break 
        end                
        xforce = ((xc-(xxi(ii)/dx))*(yc-(yyi(ii)/dy))*Fx(xf,yf) + (-1*xf + (xxi(ii)/dx))*(-1*yf + (yyi(ii)/dy))*Fx(xc,yc) + (xc-(xxi(ii)/dx))*(-1*yf + (yyi(ii)/dy))*Fx(xf,yc) + (-1*xf + (xxi(ii)/dx))*(yc-(yyi(ii)/dy))*Fx(xc,yf));
        yforce = ((xc-(xxi(ii)/dx))*(yc-(yyi(ii)/dy))*Fy(xf,yf) + (-1*xf + (xxi(ii)/dx))*(-1*yf + (yyi(ii)/dy))*Fy(xc,yc) + (xc-(xxi(ii)/dx))*(-1*yf + (yyi(ii)/dy))*Fy(xf,yc) + (-1*xf + (xxi(ii)/dx))*(yc-(yyi(ii)/dy))*Fy(xc,yf));
        
        vx(ii) = vx(ii) + dt*xforce/m;
        vy(ii) = vy(ii) + dt*yforce/m; 
        xxi(ii) = xxi(ii) + vx(ii)*dt;
        yyi(ii) = yyi(ii) + vy(ii)*dt;  
        
        
        if mod(jj-1,targ)==0 && jj>1 
            index = (jj-1)/targ + 2;
            xx(ii,index) = xxi(ii);
            yy(ii,index) = yyi(ii);
        end
        
        
        
        
        
        
        
        
        
        if DoStability == 1
            K = zeros(4);
            K(3,1) = 1;
            K(4,2) = 1;
            K(1,3) = -Fxx(xf,yf);
            K(1,4) = -Fxy(xf,yf);
            K(2,3) = -Fxy(xf,yf);
            K(2,4) = -Fyy(xf,yf);
        
            msub = squeeze(M(ii,jj,:,:));
            %det(msub)
            %M(ii,jj+1,:,:) = msub + K*msub*sqrt(vx(ii)^2+vy(ii)^2)*dt;
            M(ii,jj+1,:,:) = msub + K*msub*dt;
        end
        
        
        if CalculateVelocity == 1
            vvxx(ii,jj) = vx(ii);
            vvyy(ii,jj) = vy(ii);
        end
        
        

        


    end
end   

xx = transpose(xx);
yy = transpose(yy);





