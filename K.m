function KK=K(Y,X)
    sigma=1000;
    [n,~]=size(X);
    KK=zeros(n);
    for i=1:n
        for j=1:n
            KK(i,j)=Y(i)*Y(j)*exp(-norm(X(i,:)-X(j,:))^2/(2*sigma^2));
        end
    end
end
