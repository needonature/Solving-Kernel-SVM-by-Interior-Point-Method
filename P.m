function result=P(w,KK,t)
    [n,~]=size(KK);
    result=t*f(w,KK)-sum(log(-h(w,n)));
end
