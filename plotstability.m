np = 25;
% xx = transpose(xx);
% yy = transpose(yy);
bob = 4000;
nu = zeros(np,bob);

for ii = 1:np
    ii
    for jj = 2:bob
        mm = squeeze(M(ii,jj,:,:));
        %det(mm)
        %r = sqrt((xx(ii,jj)-xxi(ii))^2 + (yy(ii,jj)-yyi(ii))^2);
        
        nu(ii,jj) = log(abs(trace(mm)));
        %nu(ii,jj) = stability(mm,r);
    end
end


%contourf(nu)
figure(2)
plot(transpose(nu))