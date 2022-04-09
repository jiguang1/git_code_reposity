%%
%
%  PREFORMATTED
%  TEXT
%
clear all;clc;
% reference: https://www.bilibili.com/video/BV1sL4y1j7eZ/?spm_id_from=333.788
%%

%% 一阶非线性微分方程
% 齐次微分方程求解
% 齐次微分方程y'-y=0，初始条件为y(0) = 1

syms y(x)
ode =diff(y,x)-y == 0;
cond=[];
dsolve(ode,cond) %初始值为空，求通解
cond = y(0)==1; %初始值
dsolve(ode,cond)%特解
% 非齐次微分方程
% 非线性微分方程(y')^2-y=0，初始条件为y(0) =1

syms y(x)
ode = diff(y,x)-y == x;
cond=[];
dsolve(ode,cond)
cond = y(0)==1;
dsolve(ode,cond)
% 一阶非线性微分方程
% 非线性微分方程(y')^2-y=0，初始条件为y(0) =0

syms y(x)
ode = diff(y,x)^2-y==0;
cond = [];
dsolve(ode,cond)
cond = y(0)==1;
ys(x)=dsolve(ode,cond)
%% 二阶微分方程


% 二阶常系数线性微分方程
% y''+py'+qy = 0
%
% 特征方程：r^2+pr+q = 0
%
% 通解：
%
% 1.两个不相等的实数根：y=C1e^r1+C2e^r2
%
% 2.两个相等的实数根：y =(C1+C2x)e^r1
%
% 2.两个共轭复数根：y =e^(ax)*(C1*cosbx+C2sinbx)
%
%
% 1.微分方程y''-y=0，初始条件为y(0)=2，y'(0)=0

syms f(x)
ode = diff(y,x,2)-y ==0;
cond =[];
dsolve(ode,cond)
dy = diff(y,x); %赋初值前操作
cond =[y(0)==2,dy(0)==0];
dsolve(ode,cond)
% 2.微分方程x^2*y''+x*y'-y=0

syms f(x)
ode = x^2*diff(y,x,2)+x*diff(y,x)-y ==0;
cond =[];
ys(x) = dsolve(ode,cond)
% 验证求解是否正确
% diff(ys,x)
% diff(ys,x,2)
% y''-(1-y^2)y'+2y=0，初始条件为y(0) = 2，y'(0) = 0
% 没有解析解

syms f(x)
ode = diff(y,x,2)-(1-y^2)*diff(y,x)+2*y ==0;
cond =[];
ys(x) = dsolve(ode,cond)
dy = diff(y,x); %赋初值前操作
cond =[y(0)==2,dy(0)==0];
ys(x)=dsolve(ode,cond)
% 求隐士解析解
ys(x)=dsolve(ode,cond,'implicit',true)
%% 三阶微分方程
% % y'''-y=0，初始条件为y(0)=1，y'(0)=1，y''(0)=1


syms f(x)
ode = diff(y,x,3)-y ==0;
cond =[];
dsolve(ode,cond)
dy = diff(y,x); %赋初值前操作
d2y = diff(y,x,2);
cond =[y(0)==1,dy(0)==1,d2y(0)==1];
dsolve(ode,cond)
% ay'''-by=0，初始条件为y(0)=c，y'(0)=d，y''(0)=e

syms f(x) a b c d e
ode = a*diff(y,x,3)-b*y ==0;
cond =[];
dsolve(ode,cond)
dy = diff(y,x); %赋初值前操作
d2y = diff(y,x,2);
cond =[y(0)==c,dy(0)==d,d2y(0)==e];
ys(x) = dsolve(ode,cond)
simplify(ys(0)) %验证
%% 解微分方程组

% x(t)'=x+y
% y(t)'=-x+y
% 初始条件x(0) = 0,y(0)=1

syms x(t) y(t)
ode1 = diff(x,t)==x+y;
ode2 = diff(y,t)==-x+y;
cond1 = x(0)==0;
cond2 = y(0)==1;
ode = [ode1,ode2];
cond = [cond1,cond2];
[xs(t),ys(t)] = dsolve(ode,cond)
%% 微分方程数值解法ode45
% 1. y'-y=x，初始条件为 y(0) = 1；

tspan =[0,5];
y0 = 1;
opts= odeset('RelTol',1e-2,"AbsTol",1e-4);
[x,y] = ode45(@Ode1,tspan,y0,opts)

% 2. y''-(1-y^2)y'+2y=0，初始条件为 y(0) = 2，y'(0) = 0；
% 注：ode45只能求一阶微分方程，对于二阶，首先需要将其转换为一阶方程组，然后求解
%
% 因此，上式可转换为：
%
% _y1 = y_
%
% _y2 =y'_
%
% _*y1' = y2*_
%
% _*y2' = (1-y1^2)y2 -2y1*_

tspan =[0,5];
y0 = [2,0];
opts= odeset('RelTol',1e-2,"AbsTol",1e-4);
[x,y] = ode45(@Ode2,tspan,y0,opts)
sol = ode45(@Ode2,tspan,y0,opts)
% 微分方程数值解的结构体用于新值计算和扩展
% deval
%
% odextend

tspan =[0,5];
y0 = [2,0];
opts= odeset('RelTol',1e-2,"AbsTol",1e-4);
[x,y] = ode45(@Ode2,tspan,y0,opts);
sol = ode45(@Ode2,tspan,y0,opts);
x = 1;
y1 = deval(sol,x)
%%
function dy= Ode1(x,y)
dy = x+y;
end

function dy= Ode2(x,y)
dy = zeros(2,1);
dy(1) = y(2);
dy(2) = (1-y(1)^2)*y(2)-2*y(1);
end
%%
%
%
%