clear;clc;

%% 基本参数设置
m=1;Z1=1;       %入射粒子质量(单位u),电荷量
M=12;Z2=6;      %原子核质量(单位u),电荷量
% m=4;Z1=2;        %α
% M=58;Z2=28;      %58Ni

mu=m*M/(M+m);   %质心系折合质量
Z=Z1*Z2;        %电荷量乘积
k1=1.44;        %系数k1=e^2/(4*pi*varepsilon) 单位fm*MeV
V0=56;  %V0=115;
Rc=1.17*M^(1/3);a=1.5; %Woods-Saxon势参数
Ri=1.17*m^(1/3);              %入射粒子半径


%% 入射参数设置
E0=10;   %E0=139;    %入射能量，单位MeV
v0=sqrt(2*E0/mu);   %入射的经典速度
y0_list=-100+1i.*linspace(0,100,100000);  %设置入射初始坐标(虚部为b)
t_range=800/v0;                         %设置时间范围


%% 计算散射微分截面
db=y0_list(2)-y0_list(1);       %若非等距分布需更改
theta=zeros(size(y0_list));
for ii=1:length(y0_list)        %θ
    theta(ii)=scattering([y0_list(ii) v0],t_range,mu,k1,Z,V0,Rc,Ri,a);
end

Theta=abs(theta);       %θ全取正值
thetalist=linspace(min(Theta),max(Theta),1000);     %θ等距插值
dtheta=thetalist(2)-thetalist(1);                   %dθ
maxpoint=find(diff(sign([0,diff(Theta)]))==-2);
minpoint=find(diff(sign([0,diff(Theta)]))==+2);
extreme=unique([1,maxpoint,minpoint,length(Theta)]);%θ极值点位置

sigma=zeros(size(thetalist));
for jj=1:length(extreme)-1      %分段为单值函数
    bb=imag(y0_list(extreme(jj):extreme(jj+1)));
    THeta=Theta(extreme(jj):extreme(jj+1));
    col=find(thetalist(1,:)>=min(THeta)&thetalist(1,:)<=max(THeta));
    THETA=thetalist(1,col);             %只选取该段函数相应的θ插值点
    b=interp1(THeta,bb,THETA,"spline"); %插值
    db=[0,diff(b)];
    sigma(1,col)=sigma(1,col)+b.*abs(db)./(sin(THETA).*dtheta);  %σ
end


%% 计算卢瑟福散射微分截面
sigmaR=Rutherford(k1,Z,E0,thetalist);
ratio=sigma./sigmaR;

%% 绘图
figure(1)   %b-θ图
plot(rad2deg(theta),imag(y0_list))
ylabel('b(fm)')
xlabel('θ(°)')

figure(2)   %σ-θ图
plot(rad2deg(thetalist),sigma)
xlabel('θ(°)')
ylabel('σ(fm^2)')

figure(3)   %σ/σR-θ图
plot(rad2deg(thetalist),ratio,'DisplayName',...
    append('p+^{12}C',newline,'E=',num2str(E0),'MeV'))
xlabel('θ(°)')
ylabel('σ/σ_R')
ylim([1e-3,10])
set(gca, 'YScale', 'log')
xlim([0,max(rad2deg(thetalist))])
legend

%% 运动微分方程
function theta=scattering(y0,t_range,mu,k1,Z,V0,Rc,Ri,a)
tspan=[0 t_range];
options=odeset('RelTol',1e-6);   %设置容差范围
[ttraj,ytraj]=ode45(@(t,y)centralforce(t,y,mu,k1,Z,V0,Rc,Ri,a),...
                    tspan,y0,options);
theta2pi=angle(ytraj(end,2)); 
theta=theta2pi*(theta2pi<=pi)+(theta2pi-2*pi)*(theta2pi>pi);%转换坐标
end

function dy=centralforce(t,y,mu,k1,Z,V0,Rc,Ri,a)
dy=zeros(2,1);
dy(1)=y(2);
r=abs(y(1));
h=r/10000;
F=-diff(cenV([r-h,r+h],k1,Z,V0,Rc,Ri,a))/2/h;
dy(2)=F/mu*sign(y(1));
end

function V=cenV(r,k1,Z,V0,Rc,Ri,a)
R=Rc+Ri;
V=(k1*Z./r).*(r>=R)+(k1*Z*(3-r.^2/R^2)/2/R).*(r<R)-V0./(1+exp((r-Rc)./a));
%库伦势+Woods-Saxon势
end

%% 卢瑟福散射公式
function sigmaR=Rutherford(k1,Z,E0,theta)
a=k1*Z/E0;
sigmaR=a^2./(16*(sin(theta/2)).^4);
end
