%Mutual inductance

function Mab = MutualInd(ra,rb,h)

    mu=4*pi*1e-7;
    k_num=4*ra*rb;
    k_den = h^2+(ra+rb)^2;
    k=sqrt(k_num/k_den);

    [Ke,Ee]=ellipke(k*k);

    Mab=mu*sqrt(ra*rb)*((2/k-k)*Ke-2/k*Ee);

end



