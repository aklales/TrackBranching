function[A] = AbsorbingBoundary(M,N)

A = zeros(M,N);
B = zeros(M,N);

%number of points that should be affected at the boundary:
npad = 20;

for ii = 1:M
    for jj = 1:N

        if ii<npad
            A(ii,jj) = (ii-npad)^2/100; %half oscillator in x direction

            
        elseif ii>(M-npad)
            A(ii,jj) = ((ii-(M-npad))/10)^2;


         end
    end
end


for ii = 1:M
    for jj = 1:N

        if jj<npad
            B(ii,jj) = (jj-npad)^2/100; %half oscillator in x direction

            
        elseif jj>(N-npad)
            B(ii,jj) = (jj-(N-npad))^2/100;

         end
    end
end


A = A+B;