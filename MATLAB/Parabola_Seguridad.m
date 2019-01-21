clearvars -except Ts xs ys giro_m tout Vx_sim Vy_sim
close all
clc

% Angulo de giro
ang = 0.389;

% Dificultad
dif = 2;

% Altura
altura = 1;

% Altura del carril
ini = 0.16;
inc = 0.0457;
for j = 1:4
    alturas(j) = ini - (j-1)*inc;
end
for j = 5:8
    alturas(j) = -alturas(9-j);
end
alt = alturas(altura+1);
       
global c;
global s;
global t;
    
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
    global a;
    global b;
    global xF;
    global yF;

    % Valores representativos de la parábola
    global A;
    global B;
    global C;
    global D;
    global E;
    global F;
    Calculo_Constantes(ang, alt);

    % Punto del vértice
    x_ver = 0.16*s;
    y_ver = 0.16*c;

    % Derivada e inclinación del vértice
    der = -(2*A*x_ver+B*y_ver+D)/(B*x_ver+2*C*y_ver+E);
    inclinacion = atan(der);    % Angulo de inclinación en el vertice

    % Calculo de una posición de la parábola
    x = 0;
    y = Calcula_y(x);

    % Derivada e inclinación de un punto de la parábola
    der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
    alpha = atan(der);

    g = 9.8;
    R = 0.2325/2;
    m = 7.5e-3;

    % Estudio de la posición, velocidad y aceleración
    V = 0;
    x = -0.084;
    ace = 0;
    Vx=0;
    Vy=0;
    ace_x=0;
    temp = 0;
    j = 0;
    salir = 0;
    pose = 1;
%     ang_aux = ang;
    x_1 = x;
    while ((salir == 0)&&(abs(x)<0.25))
        
        temp = 0;
        
        % Caida de la moneda desde un carril al inmediatamente inferior
        if (((x>=caida(dif,altura,1))&&(caida(dif,altura,1)>0))||((x<=caida(dif,altura,1))&&(caida(dif,altura,1)<0)))
            
            altura = altura + 1;
            alt = alturas(altura+1);
            pose = -pose;
            
            % Últimas velocidades del carril
            Vx1 = V*cos(alpha);
            Vy1 = V*sin(alpha);
            
            % Calculo del punto de recepción
            Calculo_Constantes(ang, alt);
            valor_func = 1;
            while ((valor_func>0)&&(abs(x)<0.2))
                x = x + Vx1*temp;
                y = y + Vy1*temp - 0.5*g*temp^2;
                valor_func = A*x^2+B*x*y+C*y^2+D*x+E*y+F;
                
                temp = temp + 0.001;
                Vy1 = Vy1 - g*temp;
                
                % Valores para la representación
                j = j + 1;
                x_rep(j) = x;
                y_rep(j) = y;
                V_rep(j) = sqrt(Vx1^2+Vy1^2);
                a_rep(j) = g;
                t_rep(j) = j*0.001;
                
                % Fuera del sistema
                if abs(x)>0.25
                    salir = 1;
                end
                
            end
                
            % Derivada e inclinación de un punto de la parábola
            der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
            alpha = atan(der);
    
            % Velocidades en el punto de recepción
            Vx = Vx1;
            Vy = Vx*tan(alpha);
            V = Vx*sqrt(1+(tan(alpha))^2);
            ace = 2*g*sin(-alpha)/3;
                
            % Giro en el volante
%             ang_aux = -ang_aux;
        end
        
        % La caida al siguiente nivel se encuentra a la derecha
        if pose>0 

            while ((x<caida(dif,altura,1))&&(salir==0))

                Calculo_Constantes(ang, alt);
                
                % Nueva posición
                x = x + Vx*temp + 0.5*ace_x*temp^2;
                % Si se pasa del punto límete se interpola
                if x>caida(dif,altura,1)
                    gan = (caida(dif,altura,1)-x_1)/(x-x_1);
                    x = caida(dif,altura,1);
                else
                    gan = 1;
                end
                y = Calcula_y(x);

                % Inclinación de la nueva posición
                der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
                alpha = atan(der);

                % Nueva velocidad y sus componentes
                V = V + gan*ace*temp;
                Vx = V*cos(alpha);
                Vy = V*cos(alpha);

                % Nueva aceleración
                ace = 2*g*sin(-alpha)/3;
                ace_x = ace*cos(alpha);
                temp = temp + 0.001;

                % Valores para la representación                
                j = j + 1;
                x_rep(j) = x;
                y_rep(j) = y;
                V_rep(j) = V;
                a_rep(j) = ace;
                t_rep(j) = j*0.001;
                x_1 = x;
                    
                % Fuera del sistema
                if abs(x)>0.25
                    salir = 1;
                end

            end
            
        % La caida al siguiente nivel se encuentra a la derecha
        else
            
            while ((x>caida(dif,altura,1))&&(salir==0))

                Calculo_Constantes(ang, alt);
                
                % Nueva posición
                x = x + Vx*temp + 0.5*ace_x*temp^2;
                % Si se pasa del punto límete se interpola
                if x<caida(dif,altura,1)
                    gan = (x_1-caida(dif,altura,1))/(x_1-x);
                    x = caida(dif,altura,1);
                else
                    gan = 1;
                end
                y = Calcula_y(x);

                % Inclinación de la nueva posición
                der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
                alpha = atan(der);

                % Nueva velocidad 
                V = V + gan*ace*temp;
                Vx = V*cos(alpha);
                Vy = V*cos(alpha);

                % Nueva aceleración
                ace = 2*g*sin(-alpha)/3;
                ace_x = ace*cos(alpha);
                temp = temp + 0.001;

                % Valores para la representación
                j = j + 1;
                x_rep(j) = x;
                y_rep(j) = y;
                V_rep(j) = V;
                a_rep(j) = ace;
                t_rep(j) = j*0.001;
                x_1 = x;
                
                % Fuera del sistema
                if abs(x)>0.25
                    salir = 1;
                end

            end
        end
    end

    %% Representación gráfica
    
    % Representación modelo
    figure(1)
    axis([-0.25 0.25 -0.2 0.2])
    hold on
    % Parábolas
    for j = 1:size(alturas,2)
        Calculo_Constantes(ang, alturas(j));
        f = @(x,y) A*x.^2+B*x*y+C*y.^2+D*x+E*y+F;
        fimplicit(f, [-0.25143 0.25143 -0.25143 0.25143])
    end
    % Ejes de coordenadas
    pendiente = tan(ang);
    x = -0.25:0.01:0.25;
    y = -pendiente*x;
    plot(x,y,'--k')
    pendiente = 1/pendiente;
    y = -0.25:0.01:0.25;
    x = 1/pendiente*y;
    plot(x,y,'--k')
    % Simulación Simulink
%     plot(xs,ys,'g')
    % Puntos límite de los carriles
    for k = 1:size(caida,2)
        plot(caida(dif,k,1),caida(dif,k,2),'*')
    end
    % Histórico de posiciones
    plot(x_rep,y_rep,'or')
    title('Monza a -10º (dificultad 4)')
    xlabel('x')
    ylabel('y')
    grid on
    hold off
    
    figure(2)
    subplot(2,2,[1,3])
    hold on
    plot(x_rep,y_rep,'b')
    plot(xs,ys,'r')
    hold off
    title('Posición de la moneda')
    xlabel('x')
    ylabel('y')
    grid on

    subplot(2,2,2)
    for k = 1:size(Vx_sim,1)
        v_sim(k) = sqrt(Vx_sim(k)^2+Vy_sim(k)^2);
        t_rep_sim(k) = 0.00155*k;
    end
    hold on
    plot(t_rep_sim,v_sim,'r')
    plot(t_rep,V_rep,'b')
    hold off
    title('Velocidad de la moneda')
    xlabel('V (m/s)')
    ylabel('y')
    grid on

    subplot(2,2,4)
    plot(t_rep,a_rep)
    title('Aceleración de la moneda')
    xlabel('t')
    ylabel('a (m/s2)')
    grid on

else
    
    % Representación gráfica
    figure(1)
    hold on
    for j = 1:size(alturas,2)
        f = @(x,y) y + 0.54*x.^2 - alturas(j);
        fimplicit(f, [-0.25143 0.25143 -0.25143 0.25143])
    end
    for j = 1:size(caida,2)
        plot(caida(dif,j,1),caida(dif,j,2),'*')
    end
    x = -0.25:0.01:0.25;
    y = 0*x;
    plot(x,y,'--k')
    plot(y,x,'--k')
    title('Monza a 0º (dificultad 1)')
    xlabel('x')
    ylabel('y')
    grid on
    hold off
    
end





