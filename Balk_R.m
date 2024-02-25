clear
%Konstanter
A = 8;
a = 1/2;
w = 3;
tol = 1e-8;
H1 = 0.5;
HS = 0.5;
H2 = 2.8464405473;

%Visualiserar svängningsfunktionen för att hitta startvärde
% hold on
fplot(@(x) A*(exp(-a*x))*cos(w*x), [0 10]);
yline(H1,'-','H1 = 0.5');
yline(H2,'-','H2 = 2.8464405473');
grid on


%UPPGIFT 2a) Newtons metod del 1

t0 = 4.9; %Startvärde

slutvillkor = false;
iterationer = 0;
t_N1 = t0;
approximations_vektor = []; %Vektor som samlar approximationer för varje iteration för att sedan visa kvadratisk konvergens
konvergensvarde = 4.5007148743063;%Används för att visa att newton konvergerar kvadratiskt genom att mäta hur stort felet blir i varje iteration från referensvärdet
while slutvillkor == false
    iterationer = iterationer+1;
    y = A*(exp(-a*t_N1))*cos(w*t_N1) - H1;
    dy_dt = -A*(exp(-a*t_N1))*((a*cos(w*t_N1))+(w*sin(w*t_N1)));
    
    t_N1 = t_N1 - (y/dy_dt);

    y = A*(exp(-a*t_N1))*cos(w*t_N1);

    approximations_vektor = horzcat(approximations_vektor, t_N1);

    if abs(y-H1)<tol
        slutvillkor = true;
    end
end
disp(['Newtons metod 1, Startvärde: ',num2str(t0)])
disp(['Newtons metod 1, Iterationer: ',num2str(iterationer)])
disp(['Newtons metod 1, tH: ',num2str(t_N1)])

%Beräknar nu s-värdet för att visa att newton konvergerar kvadratiskt
error_vektor = [];
for i = 1:length(approximations_vektor)
    error_vektor = horzcat(error_vektor, abs(approximations_vektor(i) - t_N1));
end
for j = 2:length(error_vektor)-1
    s = error_vektor(j)/(error_vektor(j-1).^2);
    disp(['s = ',num2str(s)])
end


%UPPGIFT 2b) Sekantmetoden

x0 = [4 4.5]; %Startvärden x0 och x1

slutvillkor = false;
n = 0;
x = x0;
while slutvillkor == false
    n = n+1;
    fn_1 = A*(exp(-a*x(n)))*cos(w*x(n)) - HS;
    fn = A*(exp(-a*x(n+1)))*cos(w*x(n+1)) - HS;
    
    x(n+2) = x(n+1) - fn*((x(n+1) - x(n))/(fn - fn_1));

    f = A*(exp(-a*x(n+2)))*cos(w*x(n+2));
    if abs(f-HS)<tol
        slutvillkor = true;
    end
end
disp(['Sekantmetoden, Startvärde 1: ',num2str(x0(1)),' och startvärde 2: ',num2str(x0(2))])
disp(['Sekantmetoden, Iterationer: ',num2str(n)])
disp(['Sekantmetoden, tH: ',num2str(x(end))])
%Fler iterationer på sekantmetoden på grund av superlinjär konvergens
%istället för kvadratisk konvergent som newtons metod


%UPPGIFT 2c) Newtons metod del 2

t0 = 2.1; %Startvärde

slutvillkor = false;
iterationer = 0;
t_N2 = t0;
while slutvillkor == false
    iterationer = iterationer+1;
    y = A*(exp(-a*t_N2))*cos(w*t_N2) - H2;
    dy_dt = -A*(exp(-a*t_N2))*((a*cos(w*t_N2))+(w*sin(w*t_N2)));
    
    t_N2 = t_N2 - (y/dy_dt);

    y = A*(exp(-a*t_N2))*cos(w*t_N2);
    if abs(y-H2)<tol
        slutvillkor = true;
    end
end
disp(['Newtons metod 2, Startvärde: ',num2str(t0)])
disp(['Newtons metod 2, Iterationer: ',num2str(iterationer)])
disp(['Newtons metod 2, tH: ',num2str(t_N2)])


%UPPGIFT 2d)

%Då H är 0.5:
y1 = A*(exp(-a*t_N1))*cos(w*t_N1);
delta1 = abs(y1-H1);
disp(['delta för Newtons metod 1: ',num2str(delta1)])

%Då H är 2.8464405473:
y2 = A*(exp(-a*t_N2))*cos(w*t_N2);
delta2 = abs(y2-H2);
disp(['delta för Newtons metod 2: ',num2str(delta2)])