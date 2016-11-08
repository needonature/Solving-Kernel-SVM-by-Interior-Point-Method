function s=backtrack(z,delta_z,Y,n,t,KK)
    alpha=0.01;
    beta=0.5;
    w=z(1:n);
    u=z(n+1:2*n);
    v=z(2*n+1:3*n);
    lamda=z(end);
    u_=[u;v];
    delta_w=delta_z(1:n);
    delta_u=delta_z(n+1:2*n);
    delta_v=delta_z(2*n+1:3*n);
    delta_u_=[delta_u,delta_v];

    ii=find(delta_u_<0);
    s_max=min(1,min(-u_(ii)./delta_u_(ii)));
    s=s_max*0.99;
    w_plus=w+s*delta_w;
    while sum(h(w_plus,n)>=0)~=0
        s=beta*s;
        w_plus=w+s*delta_w;
    end
    z_plus=z+s*delta_z;
    while norm(r(z_plus,Y,n,t,KK))>(1-alpha*s)*norm(r(z,Y,n,t,KK))
        s=beta*s;
        z_plus=z+s*delta_z;
    end
end
