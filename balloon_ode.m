function dydt = balloon_ode(t,y,Mu)
% Based on:
% Mildner et al (2001) A qualitative test of the balloon model for BOLD-based MR signal changes at 3T. Magnetic resonance in medicine 46 (5) 891-9
% y(1): q (HbR)
% y(2): v (CBV)
% y(3): p (HbO) + q (HbR)
%
% This is the function to be called by ode45

% t1 = floor(t);
% t2 = ceil(t);
% if(t1==t2)
%     fin_t = Mu.fin(t);
% else
%     fin_t = interp1([t1 t2], [Mu.fin(t1) Mu.fin(t2)], t);
% end


tt = Mu.tt;
E0 = Mu.E0;
tao0 = Mu.tao0;
taov = Mu.taov;
alpha = Mu.alpha;
tau_s = Mu.tau_s;
tau_f = Mu.tau_f;
epsilon= Mu.epsilon;

% 
if any(tt == t)
    t1=find(t==tt);
    na = Mu.na(t1);
else
    [~,t1]=min(abs(t-tt));
    na = Mu.na(t1);
    
    %Mu.fin(t1),tt(t1),t
end
%     

q = y(1);
v = y(2);
fin_t = y(3);
s = y(4);
dydt = [0;0;0;0];

E_t = 1-(1-E0).^(1/fin_t);
dydt(1) = 1/tao0*(fin_t*E_t/E0-1/(tao0+taov)*(tao0*v.^(1/alpha)+taov*fin_t)*q/v);  %fin_t/tao0 * (E_t/E0 - q/v) + 1/taov * (fin_t - v.^(1/alpha))*q/v;
dydt(2) = 1/tao0*(fin_t-1/(tao0+taov)*(tao0*v.^(1/alpha)+taov*fin_t));%1/taov * (fin_t - v.^(1/alpha));
dydt(3) = 1*s;
dydt(4) = epsilon*na - s/tau_s - (fin_t-1)/tau_f;
%dydt(3) = fin_t/tao0 * ((1-E_t)/E0 - h/v) + 1/taov * (fin_t - v.^(1/alpha))*h/v;
