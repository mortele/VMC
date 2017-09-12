%% Clean up
close all;
clear variables;
clc;
format;


%% Parameters
N = 1000;
a = 1.843;
x = logspace(0,log10(5),N)-1;
x = x';


%% Wave function
norm1 = sqrt(a^3/pi);
norm2 = sqrt(a^3/(8*pi));
s1 = @(x) norm1 .* exp(-1.0*a.*x);
s2 = @(x) norm2 .* (1 - a * 0.5*x) .* exp(-a * 0.5 .* x);
y = s1(x);
figure(1);
plot(x,y,'k--','DisplayName','1s');
hold on;


%% Primitive function
g = @(c,a,x) c*exp(-a*x.^2);


%% STO-1G
f = fit(x,y,'c1*exp(-a1*x.^2)',...
        'StartPoint',randn([2 1]),...
        'Lower',[0 0],...
        'Upper',[1e4 1e4],...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e4,...
        'MaxFunEvals',1e4,...
        'DiffMinChange',1e-8,...
        'Robust','LAR')
c11 = f.c1;
a11 = f.a1;
G1 = @(x) g(c11,a11,x);
plot(x,G1(x),'DisplayName','STO-1G');

%% STO-2G
f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)',...
        'StartPoint',randn([4 1]),...
        'Lower',[0 0 0 0],...
        'Upper',[1e4 1e4 1e4 1e4],...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e4,...
        'MaxFunEvals',1e4,...
        'DiffMinChange',1e-8,...
        'Robust','LAR')
c12 = f.c1;
c22 = f.c2;
a12 = f.a1;
a22 = f.a2;
G2 = @(x) g(c12,a12,x) + g(c22,a22,x);
plot(x,G2(x),'DisplayName','STO-2G');

%% STO-3G
f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)',...
        'StartPoint',randn([6 1]),...
        'Lower',[0 0 0 0 0 0],...
        'Upper',[1e4 1e4 1e4 1e4 1e4 1e4],...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e4,...
        'MaxFunEvals',1e4,...
        'DiffMinChange',1e-8,...
        'Robust','LAR')
c13 = f.c1;
c23 = f.c2;
c33 = f.c3;
a13 = f.a1;
a23 = f.a2;
a33 = f.a3;
G3 = @(x) g(c13,a13,x) + g(c23,a23,x) + g(c33,a33,x);
plot(x,G3(x),'DisplayName','STO-3G');

%% STO-4G
f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)',...
        'StartPoint',randn([8 1]),...
        'Lower',[0 0 0 0 0 0 0 0],...
        'Upper',[1e4 1e4 1e4 1e4 1e4 1e4 1e4],...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e5,...
        'MaxFunEvals',1e5,...
        'DiffMinChange',1e-8,...
        'Robust','LAR')
c14 = f.c1;
c24 = f.c2;
c34 = f.c3;
c44 = f.c4;
a14 = f.a1;
a24 = f.a2;
a34 = f.a3;
a44 = f.a4;
G4 = @(x) g(c14,a14,x) + g(c24,a24,x) + g(c34,a34,x) + g(c44,a44,x);
plot(x,G4(x),'DisplayName','STO-4G');

%% STO-5G
f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)',...
        'StartPoint',randn([10 1]),...
        'Lower',[0 0 0 0 0 0 0 0 0 0],...
        'Upper',[1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4],...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e4,...
        'MaxFunEvals',1e4,...
        'DiffMinChange',1e-8,...
        'Robust','LAR')
c15 = f.c1;
c25 = f.c2;
c35 = f.c3;
c45 = f.c4;
c55 = f.c5;
a15 = f.a1;
a25 = f.a2;
a35 = f.a3;
a45 = f.a4;
a55 = f.a5;
G5 = @(x) g(c15,a15,x) + g(c25,a25,x) + g(c35,a35,x) + g(c45,a45,x) + g(c55,a55,x);
plot(x,G5(x),'DisplayName','STO-5G');

%% STO-6G
f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)',...
        'StartPoint',randn([12 1]),...
        'Lower',[0 0 0 0 0 0 0 0 0 0 0 0],...
        'Upper',[1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4],...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e4,...
        'MaxFunEvals',1e4,...
        'DiffMinChange',1e-8,...
        'Robust','LAR')
c16 = f.c1;
c26 = f.c2;
c36 = f.c3;
c46 = f.c4;
c56 = f.c5;
c66 = f.c6;
a16 = f.a1;
a26 = f.a2;
a36 = f.a3;
a46 = f.a4;
a56 = f.a5;
a66 = f.a6;
G6 = @(x) g(c16,a16,x) + g(c26,a26,x) + g(c36,a36,x) + g(c46,a46,x) + g(c56,a56,x) + g(c66,a66,x);
plot(x,G6(x),'DisplayName','STO-6G');




% %% GGGGGGGGGGGGG
% g = @(c,a,x) c*exp(-a*x.^2);
% %% G111111111111
% a11 = 8.549;
% c11 = 2.222;
% G1 = @(c,a,x) g(c11,a11,x);
% 
% 
% 
% %% G555555555555
% a1 = 0.7196;
% a2 = 38.68;
% a3 = 10.21;
% a4 = 1.584;
% a5 = 3.721;
% c1 = 0.02483;
% c2 = 0.7012;
% c3 = 0.7402;
% c4 = 0.2209;
% c5 = 0.5587;
% G5 = @(x) g(c1,a1,x)+g(c2,a2,x)+g(c3,a3,x)+g(c4,a4,x)+g(c5,a5,x);
% 
% nor = @(a) (2*a/pi).^(3./4);
% G   = @(c,a,x) nor(a) * c * exp(-a.*x.^2); 
% GG  = @(c,a,x) nor(a) * c * exp(-a.*x.^2);
% 
% g4_2 = @(x) G(0.4301280,11.5356690,x)       + G(0.6789140,2.0533430,x);
% g4_3 = @(x) G(0.15432897,30.1678710,x)      + G(0.53532814,5.4951153,x)       + G(0.44463454,1.4871927,x);
% g4_6 = @(x) G(0.00916359628,312.87049370,x) + G(0.04936149294,57.364462530,x) + G(0.16853830490,16.048509400,x) + G(0.37056279970,5.5130961190,x) + G(0.41649152980,2.1408965530,x) + G(0.13033408410,0.8817394283,x);
% 
% g1_2 = @(x) G(0.430128498,1.309756377,x)    + G(0.678913531,0.233135974,x);
% g1_3 = @(x) G(0.15432897,3.42525091,x)      + G(0.53532814,0.62391373,x)      + G(0.44463454,0.16885540,x);
% g1_6 = @(x) G(0.00916359628,35.52322122,x)  + G(0.04936149294,6.513143725,x)  + G(0.16853830490,1.822142904,x) + G(0.37056279970,0.625955266,x) + G(0.41649152980,0.243076747,x) + G(0.13033408410,0.100112428,x);
% 
% g1_szabo1 = @(x) GG(1,0.270950,x);
% g1_szabo2 = @(x) GG(0.678914,0.151623,x) + GG(0.430129,0.851819,x);
% g1_szabo3 = @(x) GG(0.444635,0.109818,x) + GG(0.535328,0.405771,x); + GG(0.154329,2.22766,x);
% 
% dx = mean(diff(x));
% Is1   = sum(dx * x.^2.*s1  (x).^2)
% Is11  = sum(dx * x.^2.*s11 (x).^2)
% Ig4_2 = sum(dx * x.^2.*g4_2(x).^2)
% Ig4_3 = sum(dx * x.^2.*g4_3(x).^2)
% Ig4_6 = sum(dx * x.^2.*g4_6(x).^2)
% 
% Ig1_2 = sum(dx * x.^2.*g1_2(x).^2)
% Ig1_3 = sum(dx * x.^2.*g1_3(x).^2)
% Ig1_6 = sum(dx * x.^2.*g1_6(x).^2)
% 
% Ig1_szabo1 = sum(dx * x.^2.*g1_szabo1(x).^2)
% Ig1_szabo2 = sum(dx * x.^2.*g1_szabo2(x).^2)
% Ig1_szabo3 = sum(dx * x.^2.*g1_szabo3(x).^2)
% 
% y2 = g1_6(x);
% 
% integral = @(f,x) f(x) .* x.^2;
% 
% plot(x,integral(s1,x),'k--', 'DisplayName', '1s')
% grid('on')
% hold('on')
% 
% plot(x,integral(g4_2,x), 'DisplayName', 'STO-2G')
% plot(x,integral(g4_3,x), 'DisplayName', 'STO-3G')
% plot(x,integral(g4_6,x), 'DisplayName', 'STO-6G')
% plot(x,integral(g1_2,x),   'DisplayName', 'STO-2G')
% plot(x,integral(g1_3,x),   'DisplayName', 'STO-3G')
% plot(x,integral(g1_6,x),   'DisplayName', 'STO-6G')
% plot(x,integral(g1_szabo1,x),  'DisplayName', 'Szabo STO-1G')
% plot(x,integral(g1_szabo2,x),  'DisplayName', 'Szabo STO-2G')
% plot(x,integral(g1_szabo3,x),  'DisplayName', 'Szabo STO-3G')
% 
% 
% %h = legend('Hydrogen','STO-3G','STO-6G');
% %h = legend('Hydrogen','STO-2G');
% %set(h, 'FontSize', 16,'interpreter','latex');
% h = legend('show');
% set(h,'FontSize', 18, 'interpreter', 'latex');
% %axis([0 8 0 5])
% 
% 
% 
% 
% 
