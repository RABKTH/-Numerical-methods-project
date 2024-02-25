clear
%Första raden är P1, andra raden är P2, osv.
A = [175 950; 410 2400; 675 1730];
B = [160 1008; 381 2500; 656 1760];
LA = [60; 75; 42];
LB = [45; 88; 57];
P0 = [0 0];
P4 = [1020 0];
tol = 1e-12;

%UPPGIFT 3a) MÅSTE FORTFARANDE VISA ATT DEN KONVERGERAR KVADRATISKT

%Ritar först ut cirklarna för att hitta lämpligt startvärde:
hold on
yline(0);
for n = 1:3 %Antal punkter som ska vad med i ekvationen är 3
    f1 = @(x, y) ((x - A(n,1)).^2) + ((y - A(n,2)).^2) - (LA(n).^2);
    f2 = @(x, y) ((x - B(n,1)).^2) + ((y - B(n,2)).^2) - (LB(n).^2);
    
    x = linspace(-100, 1100, 1000);
    y = linspace(-100, 2800, 1000);
    [X, Y] = meshgrid(x, y);
    Z1 = f1(X, Y);
    Z2 = f2(X, Y);
    
    contour(X, Y, Z1, [0 0], 'LineWidth', 2);
    contour(X, Y, Z2, [0 0], 'LineWidth', 2);
    grid on;
    xlabel('x-axeln');
    ylabel('y-axeln');
    title('Plot över cirklarna som ges i uppgiften');
end


%Använder nu Newtons metod för att approximera P1, P2, och P3:
startvarden = [204 1002; 458 2458; 712 1750]; %Fås manuellt genom att betrakta cirklarna
svarvektor = [];
for n = 1:3
    slutvillkor = false;
    iterationer = 0;
    x = startvarden(n,:);
    
    while slutvillkor == false
        iterationer = iterationer+1;
        
        %Räknar ut F:
        f1 = ((x(iterationer,1) - A(n,1)).^2) + ((x(iterationer,2) - A(n,2)).^2) - (LA(n).^2);
        f2 = ((x(iterationer,1) - B(n,1)).^2) + ((x(iterationer,2) - B(n,2)).^2) - (LB(n).^2);
        F = [f1; f2];
    
        %Räknar ut Jacobimatrisen DF:
        df1_dx = 2*(x(iterationer,1) - A(n,1));
        df1_dy = 2*(x(iterationer,2) - A(n,2));
        df2_dx = 2*(x(iterationer,1) - B(n,1));
        df2_dy = 2*(x(iterationer,2) - B(n,2));
        DF = [df1_dx df1_dy; df2_dx df2_dy];
        
    
        s = DF\(-F);
        x(iterationer+1,:) = x(iterationer,:) + transpose(s);
    
        if norm(s)<tol
            slutvillkor = true;
        end
    end
    svarvektor = vertcat(svarvektor, x(iterationer,:));
    disp(['P',num2str(n),', Newtons metod, Startvärde: ',num2str(startvarden(n,:))])
    disp(['P',num2str(n),', Newtons metod, Iterationer: ',num2str(iterationer)])
    disp(['P',num2str(n),', Newtons metod, x=',num2str(x(iterationer,1)),', y=',num2str(x(iterationer,2))])

    %Beräknar nu s-värdet för att visa att newton konvergerar kvadratiskt
    error_vektor = zeros(length(x),2);
    for i = 1:length(x)
        error_vektor(i,:) = abs(x(i,:) - x(iterationer,:));
    end
    for j = 2:length(error_vektor)-1
        s = norm(error_vektor(j,:))/(norm(error_vektor(j-1,:)).^2);
        disp(['P',num2str(n),', s = ',num2str(s)])
    end
end


%UPPGIFT 3b)

%p = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4

%Använder M*c=P_y_vektor där c är koefficientvektorn, M är polynomets monom, och P_y_vektor är y-värdena för de olika x-värdena

M = zeros(5,5);
j = 0;
P_x_vektor = [P0(1,1) svarvektor(1,1) svarvektor(2,1) svarvektor(3,1) P4(1,1)];
for i = P_x_vektor
    M(j+1,:) = [1, i, i.^2, i.^3, i.^4];
    j=j+1;
end

P_y_vektor = [P0(1,2) svarvektor(1,2) svarvektor(2,2) svarvektor(3,2) P4(1,2)];

c = M\transpose(P_y_vektor);

fplot(@(x) c(1) + c(2)*x + c(3)*x.^2 + c(4)*x.^3 + c(5)*x.^4)
for k = 1:5
    plot(P_x_vektor(k),P_y_vektor(k),'o')
end
disp(c)
