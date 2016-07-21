function [Gamma_salida] = wf(H,SNRlin)

r=rank(H);
p=1;
[U,D,V] = svd(H);
lambda = diag(D).^2;
Gamma=zeros(1,length(lambda));
index=[1:r];
index_temp=index;

while p<r
    irp = 1:r-p+1;
    temp = sum(1./lambda(index_temp(irp)));
    mu = 1/(r-p+1)*(1+1/SNRlin*temp);
    Gamma(index_temp(irp)) = mu - 1./(SNRlin*lambda(index_temp(irp)));
    if min(Gamma(index_temp))<0
        i=find(Gamma==min(Gamma));
        ii=find(index_temp==i);
        index_temp2=[index_temp([1:ii-1]) index_temp([ii+1:end])];
        clear index_temp;
        index_temp = index_temp2;
        p=p+1;
        clear Gamma;
    else
        p=r;
    end
end
Gamma_t = zeros(1,length(lambda));
Gamma_t (index_temp) = Gamma(index_temp);
Gamma_salida = Gamma_t;

end