function[output] = LocalForce(xx,yy,GaussianCenters,sigma_rand,b,sigma_QPC,sigma_gap,Ly,DerivativeSwitch)

%% This function calculates the potential, force, or curvature at a point
% For potential, DerivativeSwitch = 0; for Force, DerivativeSwitch = 1; for
% curvature, DerivativeSwitch = 2 (or anything else)
radius = 5*sigma_rand;



if DerivativeSwitch == 0
    %calculate potential
    
    if (abs(xx - 2*sigma_QPC) < 5*sigma_QPC) %if close to QPC, include QPC potential
        V_QPC = b*(exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)) - exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2)));
    else
        V_QPC = 0;
    end
    
    V_rand = 0;
    
    for ii = 1:length(GaussianCenters)
        if sqrt((xx - GaussianCenters(1,ii))^2 + (yy - GaussianCenters(2,ii))^2) < radius
            x0 = GaussianCenters(1,ii); 
            y0 = GaussianCenters(2,ii);
            A = GaussianCenters(3,ii);
            
            r = (x0-xx)^2 + (y0-yy)^2;
            g00 = A*exp(-r/(2*sigma_rand^2));
            
            V_rand = V_rand + g00;
        end
    end
    
    
    output = V_QPC + V_rand;
elseif DerivativeSwitch == 1
    %calculate force
    
    if (abs(xx - 2*sigma_QPC) < 5*sigma_QPC) %if close to QPC, include QPC potential
        Fx_QPC = b*(  (1/sigma_QPC^2)  *  (xx-2*sigma_QPC)  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2))  -  (1/sigma_QPC^2)  *  (xx-2*sigma_QPC)  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2)));
        Fy_QPC = b*(- (1/sigma_gap^2)  *  (yy-Ly/2)  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2)));
    else
        Fx_QPC = 0;
        Fy_QPC = 0;
    end
    
    Fx_rand = 0;
    Fy_rand = 0;
    
    for ii = 1:length(GaussianCenters)
        if sqrt((xx - GaussianCenters(1,ii))^2 + (yy - GaussianCenters(2,ii))^2) < radius
            x0 = GaussianCenters(1,ii); 
            y0 = GaussianCenters(2,ii);
            A = GaussianCenters(3,ii);
            
            r = (x0-xx)^2 + (y0-yy)^2;
            g00 = A*exp(-r/(2*sigma_rand^2));
            
            Fx_rand = Fx_rand +((xx-x0)/(sigma_rand^2))*g00;
            Fy_rand = Fy_rand +((yy-y0)/(sigma_rand^2))*g00;
            
        end
    end
    
    output = [(Fx_QPC + Fx_rand),(Fy_QPC + Fy_rand)];
else
    %calculate curvature
    
    if (abs(xx - 2*sigma_QPC) < 5*sigma_QPC) %if close to QPC, include QPC potential
        Fxx_QPC = b*(  (1/sigma_QPC^2)^2  *  (xx-2*sigma_QPC)^2  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2))  -  (1/sigma_QPC^2)^2  *  (xx-2*sigma_QPC)^2  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2))  - (1/sigma_QPC^2)  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2))  +  (1/sigma_QPC^2)  *   exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2)));
        Fyy_QPC = b*(- (1/sigma_gap^2)^2  *  (yy-Ly/2)^2  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2))  +  (1/sigma_gap^2)    *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2)));
        Fxy_QPC = b*(- (1/sigma_gap^2)  *  (yy-Ly/2)  *  (1/sigma_QPC^2)  *  (xx-2*sigma_QPC)  *  exp(-(xx-2*sigma_QPC)^2/(2*sigma_QPC^2)-(yy-Ly/2)^2 /(2*sigma_gap^2)));
    else
        Fxx_QPC = 0;
        Fyy_QPC = 0;
        Fxy_QPC = 0;
    end
    
    Fxx_rand = 0;
    Fyy_rand = 0;
    Fxy_rand = 0;
    
    for ii = 1:length(GaussianCenters)
        if sqrt((xx - GaussianCenters(1,ii))^2 + (yy - GaussianCenters(2,ii))^2) < radius
            x0 = GaussianCenters(1,ii); 
            y0 = GaussianCenters(2,ii);
            A = GaussianCenters(3,ii);
            
            r = (x0-xx)^2 + (y0-yy)^2;
            g00 = A*exp(-r/(2*sigma_rand^2));
            
            Fxx_rand = Fxx_rand + ((xx-x0)/(sigma_rand^2))^2*g00 - (1/(sigma_rand^2))*g00;
            Fyy_rand = Fyy_rand + ((yy-y0)/(sigma_rand^2))^2*g00 - (1/(sigma_rand^2))*g00;
            Fxy_rand = Fxy_rand + ((xx-x0)/(sigma_rand^2))*((yy-y0)/(sigma_rand^2))*g00;
            
        end
    end    
    
    
    output = [(Fxx_QPC + Fxx_rand),(Fxy_QPC + Fxy_rand),(Fyy_QPC + Fyy_rand)];
end