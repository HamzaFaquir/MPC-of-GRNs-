function [alpha_]=Steplength_Backtracking_mp(dJ_u1,Prot_mesh,nt_solve, ti,tf,u_k,u_l,alpha0,d1,J,kmax,P,P_j,m)


T=linspace(ti,tf, nt_solve);
bound=0;

mode=m;
L=1*norm(u_l-u_k)/norm(d1);

if L==Inf || isnan(L)
    L=alpha0;
end 

alpha=min(alpha0,L);

if alpha==0
    alpha=alpha0;
end 

gamma=(10^(-2));

k=0;








((dJ_u1.*d1))
while k<=kmax %&& alpha*norm(d1)>norm(aTc_l)

u_=u_l+alpha*d1;

 if u_<=bound
    u_=bound;
end

if u_>300
    u_=300;
end 
%}

u=ones(nt_solve,1)*u_;

[du1Lix, P_main]=main_equ(ti,tf,nt_solve,Prot_mesh,u,P_j);

J_ub=cost_f(P_main,P);




if J_ub<=J+alpha*gamma*((dJ_u1.*d1))
    disp('Right')
    alpha_=alpha;
    break
else 
    alpha=alpha*0.5;
    alpha_=alpha;

end 
k=k+1;
end 
end

 