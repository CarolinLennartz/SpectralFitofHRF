function [BOLD,fin]=balloon(tao0,taov,E0,fin,stimlen)

% This is to simulate using balloon model

% Refer to equation 11 in 
% Mildner et al (2001) A qualitative test of the balloon model for BOLD-based MR signal changes at 3T. Magnetic resonance in medicine 46 (5) 891-9

% parameters
SR=0.1;
alpha = 0.4;
%E0 = 0.3; %baseline O2 extraction
% RT = 0.1; %ramp time in sec
% DUR = 1; %duration in sec
% dCBF = 1; %delta CBF, change of cerebral blood flow
% stimTime = 1; %time of begin of stimulation in sec
% stimTime = round(stimTime/SR);

V0=0.03;
k1=6.7;
k2=2.73;
k3=0.57;

T = 32; %time in sec
tt=1:SR:T+1;

% if nargin<4
% f_in kernel
% fink =[[0:RT/SR]/RT*SR ones(1,DUR/SR-1) [RT/SR:-1:0]/RT*SR];%[[0:RT]/RT ones(1,DUR) [RT:-1:0]/RT]; %
% fink = fink*(dCBF);

% simulate fin
fin = ones(1,(T+1)/SR);
fin(1,3:3+stimlen)=2;

% fin = zeros(1,(T+1)/SR);%zeros(1,T/SR);
% fin(stimTime) = 1;%rand(size(stimTime));
% 
% fin = conv(fin, fink)+1;

%end
% if nargin==4
%     if min(fin)<1
%         fin=fin+1;
%     end
% end
if nargin<3
    E0=0.4;
end

%% now simulation

Mu = struct;
Mu.fin = fin;
Mu.E0 = E0;
Mu.tao0 = tao0;
Mu.taov = taov;
Mu.alpha = alpha;
Mu.tt=tt;

id=fin~=1;
change=find(id==1);

steps=[tt(1) tt(change(1)-1);tt(change(1)) tt(change(end));tt(change(end)) tt(end)];%{[tt(1) tt(change(1)-1)];[tt(change(1)) tt(change(end));[tt(change(end)) tt(end)]}%{tt(1:change(1)-1);tt(change(1):change(end));tt(change(end)+1:end)};
Y=zeros(length(tt),3);
step_t={1:change(1)-1;change(1):change(end);change(end)+1:length(tt)};
y0_new=[1; 1; 1];
%%%%%
for i=1:3
tspan = steps(i,:);%steps{i};%1: 0.1:T;%0:SR:T;%[1:T/SR];%
y0 = y0_new;
ode = @(t,y) balloon_ode(t,y,Mu);
[t,y] = ode45(ode, tspan, y0);
y0_new=squeeze(y(end,:));
% for j=1:3
Y(step_t{i},:)=interp1(t, y, tt(step_t{i}),'spline');%y;
% end
end
%%%%%%%%%%%%

%tt=1:0.1:T+1;
hbr=Y(:,1);
hbo=Y(:,2)-Y(:,1);
v=Y(:,2);

BOLD=V0.*(k1.*(1-hbr)+k2.*(1-(hbr./v))+k3.*(1-v));
% if nargin<3
%     fin = interp1([1:length(fin)], fin, tt);%fin(1:length(tt));%
% end


% q = ones(1,T/SR);
% v = ones(1,T/SR);
% 
% for t=1:T/SR
%     E(t) = 1-(1-E0).^(1/fin(t));
%     dq = fin(t)/tao0 * (E(t)/E0 - q(t)/v(t)) + 1/taov * (fin(t) - v(t).^(1/alpha))*q(t)/v(t);
%     dv =  1/taov * (fin(t) - v(t).^(1/alpha));
%     q(t+1) = q(t) + dq;
%     v(t+1) = v(t) + dv;
% end
% % 
% figure;plot(q); hold on;plot(v,'r')
% %figure;plot(t,y(:,1)); hold on;plot(t, y(:,2)-y(:,1),'r')
% BOLD=V0.*(k1.*(1-q)+k2.*(1-(q./v))+k3.*(1-v));
% tt=0:0.1:T;
% figure,plot(tt,BOLD)
% hold on,plot(0:0.1:T,fin(1:length(tt))-1)
% 

%%
% tt=1:0.1:T;
% hbr=y(:,1);
% hbo=y(:,2)-y(:,1);
% v=y(:,2);


% tt = 1:T/SR;%0:0.1:T;%[1.1:0.1:T];
% hbr = interp1(t,y(:,1),tt) - 1; %total deoxyhemoglobin q
% %hbo = interp1(t,y(:,2)-y(:,1),tt);
% hbo = interp1(t,y(:,3)-y(:,1),tt); %change of oxygenated Hemoglobin
% v = interp1(t,y(:,2),tt); %volume of the balloon

% BOLD=V0.*(k1.*(1-hbr)+k2.*(1-(hbr./v))+k3.*(1-v));

% % add noise?
% % noise = randn(size(hbr))/30 * 0;
% % noise2= randn(size(hbr))/60 * 0;
% % hbr = hbr + noise + noise2;
% % hbo = hbo - noise + noise2;
% 
% % filter
% %hbr = passfilter(hbr, [0.5 0.01], 10);
% %hbo = passfilter(hbo, [0.5 0.01], 10);
% %%
% % figure;
% % 
% % plotTraces([hbo'-0.5 hbr']*2, [1 2], stimTime*10, 'rbk');
% % c = runningCorrelation(hbr,hbo, 40);
% % hold on;
% % plot(c,'k')
% % xlim([1 300])
% 
% all plots
% figure;
% plot(tt,hbr,'b');
% hold on;
% plot(tt,hbo,'r');
% hold on;
% plot(tt,v+1,'k');
% hold on
% fin = interp1([1:length(fin)], fin, tt);%fin(1:length(tt));%
% plot(tt, fin+3, 'color', [0 0.5 0]);
% hold on
% 
% fout = 1/(1+tao0/taov)*(tao0/taov*v.^(1/alpha) + fin');
% plot(tt, fout+2, 'g');
% hold on
% e = 1-(1-E0).^(1./fin);
% plot(tt, e+1, 'color', [.7 .7 .7]);
% 
% legend({'HbR','HbO','volume','flow in','flow out', 'e'})
% %xlim([80 220])
% xlim([0 80])
% set(gcf,'color','w')
% axis off
% 
% figure,plot(BOLD)
% hold on,plot(fin-1)
