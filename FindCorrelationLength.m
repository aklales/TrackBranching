function[sigma,cf] = FindCorrelationLength(V,splitdim,dx)

dims = size(V);
mm = dims(1);
nn = dims(2);
Vsub = zeros(mm,nn);



if splitdim == 1
    Vsub(1:mm/2,:) = V(1:mm/2,:);
end

if splitdim == 2
    Vsub(:,1:nn/2) = V(:,1:nn/2);
end

cf = zeros(1,dims(splitdim)/2);


for ii = 1:dims(splitdim)/2
    if mod(ii,100)==0
        fprintf('Step %d of %d\n',ii,dims(splitdim)/2);
    end
    cf(ii) = sum(sum(Vsub.*V));
    Vsub = circshift(Vsub,splitdim);
end


HM = cf(1)/2;

ii = 1;

while cf(ii) > HM
    ii = ii + 1;
end

sigma = ii*dx/sqrt(2*log(2));