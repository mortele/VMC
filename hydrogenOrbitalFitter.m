%% Clean up
close all;
clear variables;
clc;
format;

%% Parameters
n = 6;
N = 1000;
Z = 4;
a = 3.983;
x1 = logspace(0,log10(2+1),N)-1;
x1 = x1';
x2 = logspace(0,log10(6+1),N)-1;
x2 = x2';
A  = logical([ones(1,n) zeros(1,n)]);
C  = logical([zeros(1,n) ones(1,n)]);

% Limits on the coefficients.
lowerC = -1e2;
upperC =  1e2;
lowerA =    0;
upperA =  1e3;
lower = @(n) [ones(1,n)*lowerA ones(1,n)*lowerC];
upper = @(n) [ones(1,n)*upperA ones(1,n)*upperC];


%% Wave function
% Node-less Slater type orbitals
norm1S = sqrt(a^3/pi);
norm2S = (1/4.)*sqrt(a^5/(6*pi));
s1S = @(x) norm1S .* exp(-a.*x);
s2S = @(x) norm2S .* x .* exp(-a * 0.5 .* x);
y1S = s1S(x1);
y2S = s2S(x2);

% Hydrogenic orbitals
norm1H = sqrt(a^3/pi);
norm2H = sqrt(a^3/(8*pi));
s1H = @(x) norm1H .* exp(-a.*x);
s2H = @(x) norm2H .* (1-a*x/2) .* exp(-a*x/2);
y1H = s1H(x1);
y2H = s2H(x2);

% 
%y1 = y1H;
%y2 = y2H;

y1 = y1S;
y2 = y2S;

% Prepare figures for plotting. 
figure(1);
plot(x1,y1,'k--','DisplayName','1s');
hold on;
figure(2);
plot(x2,y2,'k--','DisplayName','2s');
hold on;


%% Primitive function
g = @(c,a,x) c*exp(-a*x.^2);

% Gf1 = @(c,x) g(c(1),c(4),x) + g(c(1),c(5),x) + g(c(3),c(6),x);
% ccc = rand(6,1);
% cccc = lsqcurvefit(Gf1,ccc,x2,y2,lower(3),upper(3))
% figure(2);
% plot(x2,Gf1(cccc,x2),'DisplayName',strcat('STO-',num2str(n),'G'));

%% STO-nG
% 1s orbital 
functionalForm = {'c1*exp(-a1*x.^2)',
                  'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)',
                  'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)',
                  'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)',
                  'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)',
                  'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)',
                  'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)+c7*exp(-a7*x.^2)'};
%%
f = fit(x1,y1,functionalForm{n},...
        'StartPoint',randn([2*n 1]),...
        'Lower',lower(n),...
        'Upper',upper(n),...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e5,...
        'MaxFunEvals',1e5,...
        'DiffMinChange',1e-8,...
        'DiffMaxChange',0.1,...
        'Robust','LAR',...
        'TolFun',1e-6,...
        'TolX',1e-6)
ca1 = coeffvalues(f);
a1  = ca1(A);
c1  = ca1(C);
STO1= g(c1(1),a1(1),x1);
for i=2:n
    STO1 = STO1 + g(c1(i),a1(i),x1);
end
figure(1);
plot(x1,STO1,'DisplayName',strcat('STO-',num2str(n),'G'));
%%
% 2s orbital 
f = fit(x2,y2,functionalForm{n},...
        'StartPoint',randn([2*n 1]),...
        'Lower',lower(n),...
        'Upper',upper(n),...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e5,...
        'MaxFunEvals',1e5,...
        'DiffMinChange',1e-8,...
        'DiffMaxChange',0.1,...
        'Robust','LAR',...
        'TolFun',1e-6,...
        'TolX',1e-6)
ca2 = coeffvalues(f);
a2  = ca2(A);
c2  = ca2(C);
STO2= g(c2(1),a2(1),x2);
for i=2:n
    STO2 = STO2 + g(c2(i),a2(i),x2);
end

figure(2);
plot(x2,STO2,'DisplayName',strcat('STO-',num2str(n),'G'));


c = [c1;c2];
a = [a1;a2];

writeBasisToFile(Z,n,c,a,2);

%% Plot esthetics
for i=1:2
    figure(i);
    h = legend('show');
    set(h,'FontSize',18,'interpreter','latex');
end



% 
% %% STO-2G
% f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)',...
%         'StartPoint',randn([4 1]),...
%         'Lower',lower(2),...
%         'Upper',upper(2),...
%         'Algorithm','Trust-Region',...
%         'Normalize','off',...
%         'MaxIter',1e4,...
%         'MaxFunEvals',1e4,...
%         'DiffMinChange',1e-8,...
%         'Robust','LAR')
% c12 = f.c1;
% c22 = f.c2;
% a12 = f.a1;
% a22 = f.a2;
% c2 = [c12 c22];
% a2 = [a12 a22];
% G2 = @(x) g(c12,a12,x) + g(c22,a22,x);
% writeBasisToFile(Z,2,c2,a2);
% plot(x,G2(x),'DisplayName','STO-2G');
% 
% %% STO-3G
% f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)',...
%         'StartPoint',randn([6 1]),...
%         'Lower',lower(3),...
%         'Upper',upper(3),...
%         'Algorithm','Trust-Region',...
%         'Normalize','off',...
%         'MaxIter',1e4,...
%         'MaxFunEvals',1e4,...
%         'DiffMinChange',1e-8,...
%         'Robust','LAR')
% c13 = f.c1;
% c23 = f.c2;
% c33 = f.c3;
% a13 = f.a1;
% a23 = f.a2;
% a33 = f.a3;
% c3 = [c13 c23 c33];
% a3 = [a13 a23 a33];
% G3 = @(x) g(c13,a13,x) + g(c23,a23,x) + g(c33,a33,x);
% writeBasisToFile(Z,3,c3,a3);
% plot(x,G3(x),'DisplayName','STO-3G');
% 
% %% STO-4G
% f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)',...
%         'StartPoint',randn([8 1]),...
%         'Lower',lower(4),...
%         'Upper',upper(4),...
%         'Algorithm','Trust-Region',...
%         'Normalize','off',...
%         'MaxIter',1e5,...
%         'MaxFunEvals',1e5,...
%         'DiffMinChange',1e-8,...
%         'Robust','LAR')
% c14 = f.c1;
% c24 = f.c2;
% c34 = f.c3;
% c44 = f.c4;
% a14 = f.a1;
% a24 = f.a2;
% a34 = f.a3;
% a44 = f.a4;
% c4 = [c14 c24 c34 c44];
% a4 = [a14 a24 a34 a44];
% G4 = @(x) g(c14,a14,x) + g(c24,a24,x) + g(c34,a34,x) + g(c44,a44,x);
% writeBasisToFile(Z,4,c4,a4);
% plot(x,G4(x),'DisplayName','STO-4G');
% 
% %% STO-5G
% f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)',...
%         'StartPoint',randn([10 1]),...
%         'Lower',lower(5),...
%         'Upper',upper(5),...
%         'Algorithm','Trust-Region',...
%         'Normalize','off',...
%         'MaxIter',1e4,...
%         'MaxFunEvals',1e4,...
%         'DiffMinChange',1e-8,...
%         'Robust','LAR')
% c15 = f.c1;
% c25 = f.c2;
% c35 = f.c3;
% c45 = f.c4;
% c55 = f.c5;
% a15 = f.a1;
% a25 = f.a2;
% a35 = f.a3;
% a45 = f.a4;
% a55 = f.a5;
% c5 = [c15 c25 c35 c45 c55];
% a5 = [a15 a25 a35 a45 a55];
% G5 = @(x) g(c15,a15,x) + g(c25,a25,x) + g(c35,a35,x) + g(c45,a45,x) + g(c55,a55,x);
% writeBasisToFile(Z,5,c5,a5);
% plot(x,G5(x),'DisplayName','STO-5G');
% 
% %% STO-6G
% f = fit(x,y,'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)',...
%         'StartPoint',randn([12 1]),...
%         'Lower',lower(6),...
%         'Upper',upper(6),...
%         'Algorithm','Trust-Region',...
%         'Normalize','off',...
%         'MaxIter',1e4,...
%         'MaxFunEvals',1e4,...
%         'DiffMinChange',1e-8,...
%         'Robust','LAR')
% c16 = f.c1;
% c26 = f.c2;
% c36 = f.c3;
% c46 = f.c4;
% c56 = f.c5;
% c66 = f.c6;
% a16 = f.a1;
% a26 = f.a2;
% a36 = f.a3;
% a46 = f.a4;
% a56 = f.a5;
% a66 = f.a6;
% c6 = [c16 c26 c36 c46 c56 c66];
% a6 = [a16 a26 a36 a46 a56 a66]; 
% G6 = @(x) g(c16,a16,x) + g(c26,a26,x) + g(c36,a36,x) + g(c46,a46,x) + g(c56,a56,x) + g(c66,a66,x);
% writeBasisToFile(Z,6,c6,a6);
% plot(x,G6(x),'DisplayName','STO-6G');



