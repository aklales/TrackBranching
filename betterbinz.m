
xmax = 204.8;
ymax = 51.2;

nx_bins = 1024;
ny_bins = 256;
%binz2 = zeros(nx_bins,ny_bins);


for jj = 1:10000
    for ii = 1:10000
        %which bin is this point in?
        
        if ((xx(ii,jj) > 0) && (yy(ii,jj)>0))
            m = ceil(nx_bins*xx(ii,jj)/xmax);
        
            n = ceil(ny_bins*yy(ii,jj)/ymax);
        
        
            binz2(m,n) = binz2(m,n) + 1;
        end
    end
end
