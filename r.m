function [r_dual,r_cent,r_prim]=r(z,Y,n,t,K)
    w=z(1:n);
    u=z(n+1:2*n);
    v=z(2*n+1:3*n);
    lamda=z(end);
    r_dual=K*w-ones(n,1)+v-u+lamda*Y;
    r_cent=diag([u;v])*h(w,n)+1/t*ones(2*n,1);
    r_prim=Y'*w;
end
