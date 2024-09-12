%Cost function calculation 
function J=cost_f(P_main,P)
J=0.5*norm((P_main{end}-P))^2;
end 