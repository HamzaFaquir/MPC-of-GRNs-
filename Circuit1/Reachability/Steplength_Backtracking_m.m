function [alpha_]=Steplength_Backtracking_m(params,dJ_u1,u_k,u_l,d1,J)

%extract target 
P=params.opt.P_target;
T=linspace(params.time.Time_mesh(1),params.time.Time_mesh(2), params.time.Time_mesh(3));


x=cell(1,1);
Prot_mesh=params.protein.Prot_mesh;
x{1} = linspace(Prot_mesh(1,1),Prot_mesh(1,2),Prot_mesh(1,3) + 1);



%Step initialization
L=trapz(T,(u_l-u_k).^2)/abs(max(d1));
if L==Inf || isnan(L)
    L=0;
end

alpha=min(params.opt.alpha0,L);
if alpha==0
    alpha=params.opt.alpha0;
end 

gamma=params.opt.step_search.gamma;
beta=params.opt.step_search.beta;
kmax=params.opt.alpha_kmax;
alpha_=alpha;
upper_bound=params.opt.Inducer_upbo;
lower_bound=params.opt.Inducer_lobo;
k=0;



while k<=kmax && trapz(T,dJ_u1.*d1)<0 && alpha>0

u=u_l+alpha.*d1;

u(u>upper_bound)=upper_bound;
u(u<lower_bound)=lower_bound;


[du1Lix, P_main]=main_equ(params,u);

J_ub=cost_f(P_main,P,T,x);


if J_ub<=J+alpha*gamma*(trapz(T,dJ_u1.*d1))
    alpha_=alpha;
    break
else 
    alpha=alpha*beta;
    alpha_=alpha;
 
end 
k=k+1;
end 
end

 