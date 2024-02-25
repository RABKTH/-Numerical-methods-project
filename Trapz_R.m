clear
%UPPGIFT 4a)

%3 värden på x, N varieras ger:
hold on
b_lista = [0.11 0.32 1.14]; %b är x i grafen
N_intervall = [20 500];
a = 0;
Xi = sqrt(3/2);

for b = b_lista
    C_vektor = [];
    errorvektor = []; %Vektor som innehåller alla errorvärden för varje b-värde
    for N = N_intervall(1):N_intervall(2)
        h = (b - a)/N;
        x = zeros(1,N+1);
        %Skapar nu x-vektor som innehåller varje steg i intervallet
        for k = 1:(N+1)
            x(k) = a + ((k-1)*h);
        end
        %Summavektor summan av alla värden som ger approximationen om man
        %summerar dem
        summavektor = zeros(1,N+1);
        for n = 1:N+1
            fxn = (2/sqrt(pi))*exp(-((x(n)).^2));
            if n == 1
                summavektor(n) = (h/2)*fxn;
            elseif n==N+1
                summavektor(n) = (h/2)*fxn;
            else
                summavektor(n) = h*fxn;
            end
        end
        %Lägger nu till varje errorvärde för varje N i errorvektorn
        approximation = sum(summavektor);
        error = abs(erf(b) - approximation);
        errorvektor = horzcat(errorvektor, error);

        %Plottar nu C för motsvarande b-värde
        C = 1/(3*(pi.^0.5));
        C_vektor = horzcat(C_vektor, C*((b.^3)/(N.^2)));
        % C = ((h.^2)*(N.^2)*(8/(sqrt(pi)*exp(3/2))))/(12*(b.^2));
        % C_vektor = horzcat(C_vektor, C);
    end
    x_axeln = linspace(N_intervall(1),N_intervall(2),length(errorvektor));
    
    plot(x_axeln, errorvektor);
    plot(x_axeln, C_vektor,'--')
end

%legend(num2str(0.11),num2str(0.32),num2str(1.14))
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on


clear
%UPPGIFT 4b)

%3 värden på N, x varieras ger
figure(2)
hold on
b_intervall = [0.05 6];
N_lista = [50 120 400];
a = 0;
Xi = sqrt(3/2);

for N = N_lista
    C_vektor = [];
    errorvektor = []; %Vektor som innehåller alla errorvärden för varje b-värde
    for b = linspace(b_intervall(1),b_intervall(2),60) %600 är antalet punkter på x-axeln
        h = (b - a)/N;
        x = zeros(1,N+1);
        %Skapar nu x-vektor som innehåller varje steg i intervallet
        for k = 1:(N+1)
            x(k) = a + ((k-1)*h);
        end
        %Summavektor summan av alla värden som ger approximationen om man
        %summerar dem
        summavektor = zeros(1,N+1);
        for n = 1:N+1
            fxn = (2/sqrt(pi))*exp(-((x(n)).^2));
            if n == 1
                summavektor(n) = (h/2)*fxn;
            elseif n==N+1
                summavektor(n) = (h/2)*fxn;
            else
                summavektor(n) = h*fxn;
            end
        end
        %Lägger nu till varje errorvärde för varje N i errorvektorn
        approximation = sum(summavektor);
        error = abs(erf(b) - approximation);
        errorvektor = horzcat(errorvektor, error);
    end
    x_axeln = linspace(b_intervall(1),b_intervall(2),length(errorvektor));
    plot(x_axeln, errorvektor);
end
grid on
legend('N=50','N=120','N=400')
set(gca, 'YScale', 'log');