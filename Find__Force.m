function[Fx,Fy,Fxx,Fyy,Fxy,Pot] = Find__Force(x,y,GaussianCenters,sig)
%%  Finds the potential, the force, and it's derivatives at the point x,y due
%   to a number of Gaussians scattered over the sample



nGauss = length(GaussianCenters(1,:));


Fx = 0;
Fy = 0;
Fxx = 0;
Fyy = 0;
Fxy = 0;
Pot = 0;

%sum over the gaussians and evaluate at x,y
for ii = 1:nGauss
            
            x0 = GaussianCenters(1,ii); % actual location of gaussians
            y0 = GaussianCenters(2,ii);
            A = GaussianCenters(3,ii);
            
            r = (x0-x)^2 + (y0-y)^2;
            
            
        
            g00 = A*exp(-r/(2*sig^2));
            Pot = Pot + g00;
            Fx = Fx +((x-x0)/(sig^2))*g00;
            Fy = Fy +((y-y0)/(sig^2))*g00;
            Fxx = Fxx + ((x-x0)/(sig^2))^2*g00 - (1/(sig^2))*g00;
            Fyy = Fyy + ((y-y0)/(sig^2))^2*g00 - (1/(sig^2))*g00;
            Fxy = Fxy + ((x-x0)/(sig^2))*((y-y0)/(sig^2))*g00;
       
end    
       
