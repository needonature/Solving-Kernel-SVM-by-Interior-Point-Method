function DD=D_r(z,Y,n,K)
    w=z(1:n);
    u=z(n+1:2*n);
    v=z(2*n+1:3*n);
    lamda=z(end);
    DD=[K,-eye(n),eye(n),Y;[diag(-u);diag(v)],diag(h(w,n)),zeros(2*n,1);Y',zeros(1,2*n),0];
end
