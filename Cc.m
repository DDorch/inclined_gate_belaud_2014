function z=Cc(a,s,angle)
%-------------------------------------------------------------------------------------------------------
% Approximate solution for contraction coefficient  
%------------------------------------
z=1-(1.832e-2*angle-8.01e-5*angle^2)*(1-Ccv(a,s));
% A reprendre avec nouvel ajustement: Cc(free, angle)= (1.6632E-05*(90-angle)^2 + 2.2521E-03*(90-angle) + 0.611)

function z=Ccv(a,s)
% Contraction coefficient function as a combination of submerged and free flow, vertical gate - see Belaud et al. (2012)
Cc1=ccf(a);
% Calcul de p =proportion entre Cc en denoye et Cc en noye
if(s<(a*Cc1)) 
    p=1;
elseif(s>a) 
    p=0; 
else
    t=(1-s/a)/(1-Cc1); 
    p=-0.6212*t^3+1.7349*t^2-0.1133*t;
end
z=p*ccf(a)+(1-p)*ccs(a);

function z=ccf(a)
% Cc in free flow, vertical gate:
Cc0=0.618;
z=Cc0+0.0254*a^3+0.0261*a^2-0.0598*a;

%-------------------------------------------------------------------
function z=ccs(a)
% Cc in submerged flow, vertical gate
P=2*(0.194*a^2-0.499*a+0.308)+a;
z=(1-sqrt(1-P))/P;



