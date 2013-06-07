function[binz] = betterbinz(nx_bins,ny_bins,xx,yy,Lx,Ly,DontReload)


%   You choose nx_bins and ny_bins
%   xx,yy,Lx,Ly,DontReload all come from Propagate_interpolate
%   binz is a 2D histogram of the density of classical particles


if DontReload == 0
    binz = zeros(nx_bins,ny_bins);
end


bob = size(xx);

for jj = 1:bob(2)
    for ii = 1:bob(1)
        %which bin is this point in?
        
        if ((xx(ii,jj) > 0) && (xx(ii,jj) < Lx) && (yy(ii,jj)>0) && (yy(ii,jj)<Ly))
            
            m = ceil(nx_bins*xx(ii,jj)/Lx);
            n = ceil(ny_bins*yy(ii,jj)/Ly);
        
        
            binz(m,n) = binz(m,n) + 1;
        end
    end
end
