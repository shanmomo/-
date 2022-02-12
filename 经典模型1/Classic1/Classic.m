m=1;Z1=1;               %入射粒子质量(单位u),电荷量
M=12;Z2=6;               %原子核质量(单位u),电荷量
Z=Z1*Z2;                %电荷量乘积
k1=1.44;                %系数k1=e/(4*pi*varepsilon) 单位fm*MeV
V0=56;r0=1.17;a=0.75;   %Woods-Saxon势参数
R1=r0*(m^(1/3)+M^(1/3));
Rc=r0*M^(1/3);              %入射粒子半径

E0=1;               %入射能量，单位MeV
mu=m*M/(M+m);       %折合质量
v0=sqrt(2*E0/mu);   %入射的经典速度

y0_list=-50+1i.*linspace(2.34,2.35,5); %设置入射初始坐标
t_range=300/v0;                         %设置时间范围
box_l=15;                               %画图范围，单位fm
clines=cool(length(y0_list));           %颜色

for ii=1:length(y0_list)
y0=[y0_list(ii) v0];
tspan=[0 t_range];
options=odeset('RelTol',1e-6);
[ttraj,ytraj]=ode45(@(t,y)centralforce(t,y,mu,k1,Z,V0,unitlessR,r0,a),tspan,y0,options);
x=real(ytraj(:,1));
y=imag(ytraj(:,1));

%绘图
plot(x,y,"color",clines(ii,:),'DisplayName',['b=',num2str(imag(y0(1)))])
hold on
axis equal
ylim(box_l*[-1,1])
xlim(box_l*[-1,1])
end
scatter(0,0,"displayname","^{12}C")
set(gca,'color','#282C34')
legend

% %共振时位置与速度的关系
% figure(2)
% y01=[-50+2.3475i v0];
% tspan=[0 t_range];
% options=odeset('RelTol',1e-6);
% [ttraj1,ytraj1]=ode45(@(t,y)centralforce(t,y,mu,k1,Z,V0,unitlessR,r0,a),tspan,y01,options);
% yyaxis left
% plot(ttraj1,abs(ytraj1(:,1)))
% xlabel('t')
% ylabel('r(fm)')
% yyaxis right
% plot(ttraj1,abs(ytraj1(:,2)))
% ylabel('v(fm/s)')
% xlim([20,160])

function dy=centralforce(t,y,mu,k1,Z,V0,R1,Rc,a)
dy=zeros(2,1);
dy(1)=y(2);
r=abs(y(1));
h=r/10000;
F=-diff(cenV([r-h,r+h],k1,Z,V0,R1,Rc,a))/2/h;
dy(2)=F/mu*sign(y(1));
end

function V=cenV(r,k1,Z,V0,R1,Rc,a)
V=(k1*Z./r).*(r>=R1)+(k1*Z*(3-r.^2/R1^2)/2/R1).*(r<R1)-V0./(1+exp((r-Rc)./a));
%库伦势+Woods-Saxon势
end
