clearvars -except Ts xs ys giro_m tout
close all
clc

% Angulo de giro
ang = 0*pi/12;

% Dificultad
dif = 3;

% Altura del carril
alturas = [0.16 0.1143 0.0686 0.03 -0.03 -0.0686 -0.1143 -0.16];
alt = alturas(2);
       
% Trigonometría
c = cos(ang);
s = sin(ang);
t = tan(ang);
t = 1/t;

% Vectores de la posición de caida para distintas dificultades
caida(1,:,:) = [0 0.11429;
           0 0.06857;
           0 0.02286;
           0 -0.02286;
           0 -0.06857;
           0 -0.11429;
           -0.02554 -0.16035];
       
caida(2,:,:) = [0.06211 0.11223;
           -0.04758 0.06799;
           0.04969 0.02154;
           -0.04847 -0.02411;
           0.04406 -0.06961;
           -0.05409 -0.11585;
           -0.02554 -0.16035];
       
caida(3,:,:) = [0.12394 0.10605;
           -0.09503 0.06374;
           0.09924 0.01759;
           -0.0968 -0.02787;
           0.08803 -0.07271;
           -0.108 -0.12053;
           -0.02554 -0.16035];
       
caida(4,:,:) = [0.12394 0.10605;
           -0.14224 0.05771;
           0.14851 0.01102;
           -0.14488 -0.03412;
           0.1318 -0.07789;
           -0.108 -0.12053;
           -0.02554 -0.16035];

if ang~=0
    % Cambio a sistema de referencia no variable
    for j = 1:size(caida,1)
        for k = 1:size(caida,2)
            x = caida(j,k,1);
            y = caida(j,k,2);
            dist = sqrt(x^2+y^2);
            if y>=0
                gamma = atan(x/y);
            elseif x>=0
                gamma = pi + atan(x/y);
            else
                gamma = -pi + atan(x/y);
            end

            caida(j,k,1) = dist*sin(gamma+ang);
            caida(j,k,2) = dist*cos(gamma+ang);
        end
    end

    % Reducción de complejidad
    a = t - 1/(s*c);
    b = (alt + 0.463)/c;

    % Punto focal
    xF = (alt-0.463)*s;
    yF = (alt-0.463)*c;

    % Valores representativos de la parábola
    global A;
    global B;
    global C;
    global D;
    global E;
    global F;
    A = 1;
    B = 2*a;
    C = a^2;
    D = -2*(xF*(a^2+1)+b*a);
    E = 2*(b-(a^2+1)*yF);
    F = (a^2+1)*xF^2 + (a^2+1)*yF^2 - b^2;

    % Punto del vértice
    x_ver = 0.16*s;
    y_ver = 0.16*c;

    % Derivada e inclinación del vértice
    der = -(2*A*x_ver+B*y_ver+D)/(B*x_ver+2*C*y_ver+E);
    inclinacion = atan(-der);    % Angulo de inclinación en el vertice

    % Calculo de una posición de la parábola
    x = 0;
    y = Calcula_y(x);

    % Derivada e inclinación de un punto de la parábola
    der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
    alpha = atan(-der);

    g = 9.8;
    R = 0.2325/2;
    m = 7.5e-3;

    % Estudio de la posición, velocidad y aceleración
    V = 0;
    x = 0;
    ace = 0;
    t = 0;
    j = 0;
    while abs(x)<0.21

        % Nueva posición
        x = x + V*t + 0.5*ace*t^2;
        y = Calcula_y(x);

        % Inclinación de la nueva posición
        der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
        alpha = atan(-der);

        % Nueva velocidad y sus componentes
        V = V + ace*t;
        Vx1 = V*cos(alpha);
        Vy1 = V*sin(alpha);

        % Nueva aceleración
        ace = 2*g*sin(alpha)/3;
        t = t + 0.01;
        j = j + 1;

        % Valores para la representación
        x_rep(j) = x;
        y_rep(j) = y;
        V_rep(j) = V;
        a_rep(j) = ace;
        t_rep(j) = t;

    end

    % Representación gráfica
    figure(1)
    hold on
    for j = 1:size(alturas,2)
        f = @(x,y) A*x.^2+B*x*y+C*y.^2+D*x+E*y+F;
        fimplicit(f, [-0.25143 0.25143 -0.25143 0.25143])
    end
    plot(xs,ys)
    plot(caida(3,1,1),caida(3,1,1),'*')
    title('Parabola a 30º')
    xlabel('x')
    ylabel('y')
    grid on
    hold off
    
    figure(2)
    subplot(2,2,[1,3])
    plot(x_rep,y_rep)
    title('Posición de la moneda')
    xlabel('x')
    ylabel('y')
    grid on

    subplot(2,2,2)
    plot(t_rep,V_rep)
    title('Velocidad de la moneda')
    xlabel('t')
    ylabel('y')
    grid on

    subplot(2,2,4)
    plot(t_rep,a_rep)
    title('Aceleración de la moneda')
    xlabel('t')
    ylabel('y')
    grid on

else
    
    % Representación gráfica
    figure(1)
    hold on
    for j = 1:size(alturas,2)
        f = @(x,y) y + 0.54*x.^2 - alturas(j);
        fimplicit(f, [-0.25143 0.25143 -0.25143 0.25143])
    end
    plot(xs,ys)
    plot(caida(3,1,1),caida(3,1,1),'*')
    title('Parabola a 30º')
    xlabel('x')
    ylabel('y')
    grid on
    hold off
    
end





