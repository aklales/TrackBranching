function[A] = Z_AbsorbingBoundary(M,N,w)

A = zeros(M,N);
B = zeros(M,N);

p = 3.0;
%number of points that should be affected at the boundary:
npad = w*min(M,N);

for ii = 1:M
    for jj = 1:N

        if ii<npad
            A(ii,jj) = (ii-npad)^2/p; %half oscillator in x direction

            
        elseif ii>(M-npad)
            A(ii,jj) = (ii-(M-npad))^2/p;


         end
    end
end


for ii = 1:M
    for jj = 1:N

        if jj<npad
            B(ii,jj) = (jj-npad)^2/p; %half oscillator in x direction

            
        elseif jj>(N-npad)
            B(ii,jj) = (jj-(N-npad))^2/p;

         end
    end
end


A = A+B;