function result=g(w,KK,t)
    [n,~]=size(KK);
    C=1000;
    result=t*KK*w-t*ones(n,1)-1./w+1./(C-w);
end
