function J=cost_f(P_main,P,T,x)

%{
[r,o]=size(T);

I=zeros(1,o);

for i=1:o
    I(i)=norm((P_main{i}-P))^2;

end 
%}

J=0.5*norm(P_main{end}-P)^2;%0.5*trapz(x{1},(P_main{end}-P).^2)^2;%+0.5*trapz(T,I);

end 