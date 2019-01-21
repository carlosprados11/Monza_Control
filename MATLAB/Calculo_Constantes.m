function Calculo_Constantes(ang, alt)

    global c;
    global s;
    global t;
    global a;
    global b;
    global xF;
    global yF;
    
    % Trigonometría
    c = cos(ang);
    s = sin(ang);
    t = tan(ang);
    t = 1/t;

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

end

