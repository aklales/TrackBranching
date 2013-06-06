%%  Propagate the particles

%clear all;
D2 = 1;
m = 1;
DoStability = 0;
PEHO = 8;

np = 10000;


vel = 4;                % initial velocity of particle
time = 10000;
dt = 0.01;

%load('Potential_4096_1024_i.mat')%loads Force and coords of entries in force matrix


%nW = length(xx);
%nL = length(yy);

%xV = xx;
%yV = yy;

%Fx = Fx_rand + Fx_QPC;
%Fy = Fy_rand + Fy_QPC;

%Vrand = V;

%V = Vrand + V_QPC;



        % number of grid points along


%dx = abs(xV(2)-xV(1));  %spatial difference between adjacent points of potential
%dy = abs(yV(2)-yV(1));
%Lx = xV(end);           %real dimensions of rectangle 
%Ly = yV(end);


%[Fx,Fy] = forcify(V,dx,dy);
if DoStability ==1
%    [Fxx,Fxy] = forcify(Fx,dx,dy);
%    [Fyx,Fyy] = forcify(Fy,dx,dy);
    Fxx = Fxx_rand + Fxx_QPC;
    Fyy = Fyy_rand + Fyy_QPC;
    Fxy = Fxy_rand + Fxy_QPC;
end
%clear('Fx_rand','Fy_rand','Fxy_rand','Fyy_rand','Fxx_rand','Fx_QPC','Fy_QPC','Fxx_QPC','Fxy_QPC','Fyy_QPC','Vimag','Vrand','V_QPC');

xx = zeros(np,time); 
yy = zeros(np,time);
vvxx = zeros(np,time);
vvyy = zeros(np,time);


%%  SET UP INTIAL CONDITIONS

    

vy = zeros(np,1); %initial momenta of particles
vx =  vel*ones(np,1); %zeros(np,1); 

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


vvxx(:,1) = vx;
vvyy(:,1) = vy;
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
for ii = 1:np
    if mod(ii,10)==0
        ii
    end
    for jj = 1:time
        
        xf = floor(xx(ii,jj)/dx);
        xc = ceil(xx(ii,jj)/dx);
        yf = floor(yy(ii,jj)/dy);
        yc = ceil(yy(ii,jj)/dy);

        if ((xf < 25) || (xf > (nW-25)) || (yf < 25) || (yf > (nL-25))) 
                
            break 
        end                
        xforce = ((xc-(xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fx(xf,yf) + (-1*xf + (xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fx(xc,yc) + (xc-(xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fx(xf,yc) + (-1*xf + (xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fx(xc,yf));
        yforce = ((xc-(xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fy(xf,yf) + (-1*xf + (xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fy(xc,yc) + (xc-(xx(ii,jj)/dx))*(-1*yf + (yy(ii,jj)/dy))*Fy(xf,yc) + (-1*xf + (xx(ii,jj)/dx))*(yc-(yy(ii,jj)/dy))*Fy(xc,yf));
%         Fx = Force(xint,yint,1);%   which forces should we use?
%         if D2 == 1;
%             Fy = Force(xint,yint,2);
%         else
%             Fy = 0;
%         end
        
        vx(ii) = vx(ii) + dt*xforce/m;
        vy(ii) = vy(ii) + dt*yforce/m; 
        xx(ii,jj+1) = xx(ii,jj) + vx(ii)*dt;
        yy(ii,jj+1) = yy(ii,jj) + vy(ii)*dt;  
        
        if DoStability ==1
            K = zeros(4);
            K(3,1) = 1;
            K(4,2) = 1;
            K(1,3) = -Fxx(xf,yf);
            K(1,4) = -Fxy(xf,yf);
            K(2,3) = -Fxy(xf,yf);
            K(2,4) = -Fyy(xf,yf);
        
            M(ii,jj+1,:,:) = squeeze(M(ii,jj,:,:)) + K*squeeze(M(ii,jj,:,:))*sqrt(vx(ii)^2+vy(ii)^2)*dt;
        end
        
        

        vvxx(ii,jj) = vx(ii);
        vvyy(ii,jj) = vy(ii);

    end
end   

  
%%
% if DoStability == 1
%     save('june5.mat','xx','yy','vvxx','vvyy','M')
% else
%     save('june5.mat','xx','yy','vvxx','vvyy')
% end
% figure(5)
 %ColorSet = varycolor(np);
 %set(0,'DefaultAxesColorOrder',ColorSet);
 %plot(transpose(xx),transpose(yy),'g.','MarkerSize',1)
