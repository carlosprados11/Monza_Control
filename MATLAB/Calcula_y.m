function y = Calcula_y(x)

    global A;
    global B;
    global C;
    global D;
    global E;
    global F;
    
    y = (-(B*x+E)+sqrt((B*x+E)^2-4*C*(A*x^2+D*x+F)))/(2*C);
    if (abs(y)>0.25)
        y = (-(B*x+E)-sqrt((B*x+E)^2-4*C*(A*x^2+D*x+F)))/(2*C);
    end
        
end

