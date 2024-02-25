clear
load STHLMARLANDA.mat
w = (2.*pi)/(365.*24);

%UPPGIFT 1a)

%T = c1 + c2*t + A0*sin(w*t) + A1*cos(w*t);
%A0 = c3*cos(w*ts);
%A1 = c3*(-sin(w*ts));

%Designmatrisen, Ax=b där Td=b i vårt fall:
A = zeros(length(Td),4);
for t = 0:length(Td)-1 %Tiden går från 0 till 130813, alltså 130814 "steg"
    A(t+1,:) = [1, t, sin(w.*t), cos(w.*t)];
end

%Koefficientmatrisen, x:
AT = transpose(A);
x = (AT*A)\(AT*Td);


%UPPGIFT 1b)

%Beräknar c1, c2, c3 och ts:
c1 = x(1);
c2 = x(2);
disp(['c1 = ',num2str(c1)]);
disp(['c2 = ',num2str(c2)]);
A0 = x(3);
A1 = x(4);
c3 = sqrt((A0.^2) + (A1.^2));
disp(['c3 = ',num2str(c3)]);

ts = (1/w)*atan(-A1/A0);
disp(['ts = ',num2str(ts)]);

%Plottar:
t = (0:(length(Td)-1));
scatter(t, Td, 5)
title('Funktion och givna värden')
hold on
tinterval = [0 (length(Td)-1)];
fplot(@(t) c1 + c2*t + A0*sin(w*t) + A1*cos(w*t), tinterval)
hold off


%UPPGIFT 1c)

%Räknar ut Tmod:
tspan = (0:length(Td)-1);
Tmod = zeros(length(Td),1);
for t = tspan
    Tmod(t+1) = c1 + c2*t + A0*sin(w*t) + A1*cos(w*t);
end
Tdiff = Td-Tmod;

%Beräknar 2-normen av residualen:
residual_2_norm = norm(Tdiff,2);

%Beräknar L^2-normen:
diskret_L2_norm = norm(Tdiff,2)/(sqrt(length(Td)));


%UPPGIFT 1d)

%Normalekvation och backslash, Ax=b -> ATAx=ATb:
AT = transpose(A);
x1 = (AT*A)\(AT*Td);

%Backslash direkt:
x2 = A\Td;

%Skillnaden mellan metoderna:
x_skillnad = x1-x2;
disp(['Skillnaden i koefficienter mellan normalekvation och direkt backslash är c1=',num2str(x_skillnad(1)),', c2=',num2str(x_skillnad(2)),', A0=',num2str(x_skillnad(3)),', A1=',num2str(x_skillnad(4))]);

%Konditionstalet för matrisen för normalekvationerna:
cond1 = cond(AT*A, inf);
disp(['Konditionstalet är ',num2str(cond1)])
