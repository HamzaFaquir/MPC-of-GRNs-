function [alpha_]=Steplength_Backtracking_m(params,dJ_u1,u_k,u_l,d1,J)



iN=cell(1,1);
x=cell(1,1);
iN{1} = params.protein.Prot_mesh(1,3) + 1;
x{1} = linspace(params.protein.Prot_mesh(1,1),params.protein.Prot_mesh(1,2), iN{1});

%extract target 
P=params.opt.P_target;
T=linspace(params.time.Time_mesh(1),params.time.Time_mesh(2), params.time.Time_mesh(3));

%Intial step calculation
L=1*trapz(T,(u_l-u_k).^2)/abs(max(d1));
if L==Inf || isnan(L)
    L=0;
end

alpha=min(params.opt.alpha0,L);

if alpha==0
    alpha=params.opt.alpha0;
end 


%Extract line search parameters
gamma=params.opt.step_search.gamma;
beta=params.opt.step_search.beta;

alpha_=alpha;
k=0;


upper_bound=params.opt.Inducer_upbo;
lower_bound=params.opt.Inducer_lobo;




while k<=params.opt.alpha_kmax && trapz(T,dJ_u1.*d1)<0 && alpha>0

u=u_l+alpha.*d1;

u_b=u>upper_bound;
l_b=u<lower_bound;


u(u_b)=upper_bound;
u(l_b)=lower_bound;



[du1Lix, P_main]=main_equ(params,u);

J_ub=cost_f(P_main,P,T,x);


if J_ub<=J+alpha*gamma*(trapz(T,dJ_u1.*d1))
    disp('Right')
    alpha_=alpha;
    break
else 
    alpha=alpha*beta;
    alpha_=alpha;
 
end 
k=k+1;
end 
end

 