clearvars -except Ts xs ys giro_m tout Vx_sim Vy_sim
close all
clc

% Periodo de muestreo
Ts = 0.033;
inc_t = 0.001;
inc_t = Ts;

% Cargo el controlador
controlfuzzy = readfis('Controlador_3');

% Angulo de giro
ang = 0;

% Dificultad
dif = 4;

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
    
% Vectores de la posición de caida para distintas dificultades
caida_i(1,:,:) = [0 0.11429;
           0 0.06857;
           0 0.02286;
           0 -0.02286;
           0 -0.06857;
           0 -0.11429;
           -0.02554 -0.16035];
       
caida_i(2,:,:) = [0.06211 0.11223;
           -0.04758 0.06799;
           0.04969 0.02154;
           -0.04847 -0.02411;
           0.04406 -0.06961;
           -0.05409 -0.11585;
           -0.02554 -0.16035];
       
caida_i(3,:,:) = [0.12394 0.10605;
           -0.09503 0.06374;
           0.09924 0.01759;
           -0.0968 -0.02787;
           0.08803 -0.07271;
           -0.108 -0.12053;
           -0.02554 -0.16035];
       
caida_i(4,:,:) = [0.12394 0.10605;
           -0.14224 0.05771;
           0.14851 0.01102;
           -0.14488 -0.03412;
           0.1318 -0.07789;
           -0.108 -0.12053;
           -0.02554 -0.16035];

       
% Representación del historico
figure(1)
axis([-0.25 0.25 -0.25 0.25])
hold on
for j = 1:size(alturas,2)
    f = @(x,y) y + 0.54*x.^2 - alturas(j);
    fimplicit(f, [-0.25143 0.25143 -0.25143 0.25143])
end
for j = 1:size(caida_i,2)
    plot(caida_i(dif,j,1),caida_i(dif,j,2),'*')
end
x = -0.25:0.01:0.25;
y = 0*x;
plot(x,y,'--k')
plot(y,x,'--k')
hold off
title('Puntos por los que pasa la moneda (dificultad 4)')
xlabel('x')
ylabel('y')
grid on

% Valores iniciales
V = 0;
x = -0.11;
x_1 = x;
y = -0.54*x^2 + alt;
y_1 = y;
ang_1 = 0;
ace = 0;
Vx=0;
Vy=0;
ace_x=0;
temp = 0;
j = 0;
estado = 0;
cambio = 1;
distancia_ver = 0;
movimiento = [1 0 1 0 1 0 1];
limite = pi/6;

% Parámetros
g = 9.8;
R = 0.2325/2;
m = 7.5e-3;
umbral_ang = 1e-5;
R1 = 0.2172;
    
while (((x^2+y^2)<R1^2)&&(altura<8))
    
    caida = caida_i;
    j = j + 1;
    ang_rep(j) = ang;
    
    if ang~=0
        % Cambio a sistema de referencia no variable en el modo de
        % dificultad dado
        
        for k = 1:size(caida,2)
            x_c = caida(dif,k,1);
            y_c = caida(dif,k,2);
            dist = sqrt(x_c^2+y_c^2);
            if y_c>=0
                gamma = atan(x_c/y_c);
            elseif x_c>=0
                gamma = pi + atan(x_c/y_c);
            else
                gamma = -pi + atan(x_c/y_c);
            end

            caida(dif,k,1) = dist*sin(gamma+ang);
            caida(dif,k,2) = dist*cos(gamma+ang);
        end

        Calculo_Constantes(ang, alt);

    end
    
    % Cambio de sistema de referencia
    if estado==0
        
        dist = sqrt(x^2+y^2);
        if y>=0
            gamma = atan(x/y);
        elseif x>=0
            gamma = pi + atan(x/y);
        else
            gamma = -pi + atan(x/y);
        end

        x = dist*sin(gamma+(ang-ang_1));
        x_1=x;
    end

    % Caida de la moneda desde un carril al inmediatamente inferior
    if (estado==1)%((((x>=caida(dif,altura,1))&&(caida(dif,altura,1)>0))||((x<=caida(dif,altura,1))&&(caida(dif,altura,1)<0)))&&(estado==1))

        % Calculo del punto de recepción
        if ((ang~=0) && (altura~=7))
            Calculo_Constantes(ang, alturas(altura+2));
        end
                
        x = x + Vx*temp;
        y = y + Vy*temp - 0.5*g*temp^2;

        if ((ang==0) && (altura~=7))
            valor_func = y + 0.54*x^2 - alturas(altura+2);
        elseif altura~=7
            valor_func = A*x^2+B*x*y+C*y^2+D*x+E*y+F;
        end
        
        % Nueva velocidad en y
        Vy = Vy - g*temp;
        temp = temp + inc_t;

        % Si está en la última etapa no se calcula nada
        if altura==7
            valor_func = 1;
        end
        
        % Si toca el siguiente carril
        if valor_func <=0
            estado = 0;
            altura = altura + 1;
            alt = alturas(altura+1);
            temp = 0;
            
            % Derivada e inclinación de un punto de recepción
            if ang==0
                der = -2*0.54*x;
            else
                der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
            end
            alpha = atan(der);
        
            % Velocidades en el punto de recepción
            Vy = Vx*tan(alpha);
            V = Vx*sqrt(1+(tan(alpha))^2);
            ace = 2*g*sin(-alpha)/3;
            ace_x = ace*cos(alpha);
        end

        % Valores para la representación
        x_rep(j) = x;
        y_rep(j) = y;
        V_rep(j) = sqrt(Vx^2+Vy^2);
        a_rep(j) = g;
        t_rep(j) = j*inc_t;


    % Movimiento por un carril donde la caida se encuentra a la derecha
    elseif movimiento(altura)==1 %(((caida(dif,altura,1)>=0)&&(x<caida(dif,altura,1)))||((altura==7)&&(estado==0)))

        if ang~=0
            Calculo_Constantes(ang, alt);
        end

        % Nueva posición
        x = x + Vx*temp + 0.5*ace_x*temp^2;
        
        % Si se pasa del punto límete se interpola
        if x>=caida(dif,altura,1)
            gan = (caida(dif,altura,1)-x_1)/(x-x_1);
            x = caida(dif,altura,1);
            estado = 1;
            temp = -inc_t;
        else
            gan = 1;
        end
        
        % Calculo de y
        if ang==0
            y = -0.54*x^2+alt;
        else
            y = Calcula_y(x);
        end

        % Derivada e inclinación de un punto de la parábola
        if ang==0
            der = -2*0.54*x;
        else
            der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
        end
        alpha = atan(der);

        % Nueva velocidad y sus componentes
        V = V + gan*ace*temp;
        Vx = V*cos(alpha);
        Vy = V*sin(alpha);

        % Nueva aceleración
        ace = 2*g*sin(-alpha)/3;
        ace_x = ace*cos(alpha);
        temp = temp + inc_t;

        % Valores para la representación                
        x_rep(j) = x;
        y_rep(j) = y;
        V_rep(j) = V;
        a_rep(j) = ace;
        t_rep(j) = j*inc_t;

    % Movimiento por un carril donde la caida se encuentra a la izquierda
    else %((x>caida(dif,altura,1))&&(estado==0))

        if ang~=0
            Calculo_Constantes(ang, alt);
        end

        % Nueva posición
        x = x + Vx*temp + 0.5*ace_x*temp^2;
        % Si se pasa del punto límete se interpola
        if x<=caida(dif,altura,1)
            gan = (x_1-caida(dif,altura,1))/(x_1-x);
            x = caida(dif,altura,1);
            estado = 1;
            temp = -inc_t;
        else
            gan = 1;
        end
        
        % Calculo de y
        if ang==0
            y = -0.54*x^2+alt;
        else
            y = Calcula_y(x);
        end

        % Derivada e inclinación de un punto de la parábola
        if ang==0
            der = -2*0.54*x;
        else
            der = -(2*A*x+B*y+D)/(B*x+2*C*y+E);
        end
        alpha = atan(der);

        % Nueva velocidad 
        V = V + gan*ace*temp;
        Vx = V*cos(alpha);
        Vy = V*sin(alpha);

        % Nueva aceleración
        ace = 2*g*sin(-alpha)/3;
        ace_x = ace*cos(alpha);
        temp = temp + inc_t;

        % Valores para la representación
        x_rep(j) = x;
        y_rep(j) = y;
        V_rep(j) = V;
        a_rep(j) = ace;
        t_rep(j) = j*inc_t;

    end
    
    % Almacenar la posición anterior
    ang_1 = ang;
    x_1 = x;
    
    if estado==0 % En el movimiento por el carril (cuando cae no movemos)
       
        if movimiento(altura)==1 %((caida_i(dif,altura,1)>=0)||(altura==7))
            % Calculo de las entradas del controlador
            dist = caida(dif,altura,1) - x;
            velocidad = Vx;
            inclin = alpha;
            % Control y calculo del nuevo angulo
            ang = ang + evalfis([dist inclin velocidad],controlfuzzy);
        else
            % Calculo de las entradas del controlador
            dist = x - caida(dif,altura,1);
            velocidad = -Vx;
            inclin = -alpha;
            % Control y calculo del nuevo angulo
            ang = ang - evalfis([dist inclin velocidad],controlfuzzy);
        end
        
        if abs(ang)<umbral_ang
            ang = 0;
        end
        if ang>limite
            ang = limite;
        elseif ang<-limite
            ang = -limite;
        end
        ang

    end
    
    % Representación modelo
    figure(1)
    hold on
    plot(x,y,'or')
    hold off
    
end

%% Representación gráfica

% Representación modelo
% figure(1)
% axis([-0.25 0.25 -0.2 0.2])
% hold on
% % Parábolas
% % for j = 1:size(alturas,2)
% %     Calculo_Constantes(ang, alturas(j));
% %     f = @(x,y) A*x.^2+B*x*y+C*y.^2+D*x+E*y+F;
% %     fimplicit(f, [-0.25143 0.25143 -0.25143 0.25143])
% % end
% % % Ejes de coordenadas
% % pendiente = tan(ang);
% % x = -0.25:0.01:0.25;
% % y = -pendiente*x;
% % plot(x,y,'--k')
% % pendiente = 1/pendiente;
% % y = -0.25:0.01:0.25;
% % x = 1/pendiente*y;
% % plot(x,y,'--k')
% % % Simulación Simulink
% % plot(xs,ys,'g')
% % % Puntos límite de los carriles
% % for k = 1:size(caida,2)
% %     plot(caida(dif,k,1),caida(dif,k,2),'*')
% % end
% % Histórico de posiciones
% plot(x_rep,y_rep,'or')
% title('Puntos por los que pasa la moneda (dificultad 4)')
% xlabel('x')
% ylabel('y')
% grid on
% hold off

figure(2)
subplot(2,2,3)
hold on
plot(t_rep,ang_rep,'b')
% plot(xs,ys,'r')
hold off
title('Giro de la Monza')
xlabel('t')
ylabel('Giro (rad)')
grid on

subplot(2,2,1)
hold on
plot(x_rep,y_rep,'b')
% plot(xs,ys,'r')
hold off
title('Posición de la moneda')
xlabel('x')
ylabel('y')
grid on

subplot(2,2,2)
% for k = 1:size(Vx_sim,1)
%     v_sim(k) = sqrt(Vx_sim(k)^2+Vy_sim(k)^2);
%     t_rep_sim(k) = 0.00155*k;
% end
hold on
% plot(t_rep_sim,v_sim,'r')
plot(t_rep,V_rep,'b')
hold off
title('Velocidad de la moneda')
xlabel('t')
ylabel('V (m/s)')
grid on

subplot(2,2,4)
plot(t_rep,a_rep)
title('Aceleración de la moneda')
xlabel('t')
ylabel('a (m/s2)')
grid on



% % Representación gráfica
% figure(1)
% hold on
% for j = 1:size(alturas,2)
%     f = @(x,y) y + 0.54*x.^2 - alturas(j);
%     fimplicit(f, [-0.25143 0.25143 -0.25143 0.25143])
% end
% for j = 1:size(caida,2)
%     plot(caida(dif,j,1),caida(dif,j,2),'*')
% end
% x = -0.25:0.01:0.25;
% y = 0*x;
% plot(x,y,'--k')
% plot(y,x,'--k')
% title('Monza a 0º (dificultad 1)')
% xlabel('x')
% ylabel('y')
% grid on
% hold off


