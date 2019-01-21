function [x_choque,y_choque,x_choque_m]=solver(a,b,v,ang_pos,x_resp,y_resp,piso,timer_choque)

x_choque_m=zeros(4,1);
if (piso==1 || piso==3 || piso==5)
    x_choque=10;
else
    if (piso==2 || piso==4 || piso==6)
        x_choque=-10;
    end
end
y_choque=0;
d=x_resp;
h=y_resp;

syms xc;
if(timer_choque>=2)
     switch piso
        case 1
            pto_choque=real(double(solve((-b*xc+a*((xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h))==-0.54*(a*xc+b*(xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h).^2+0.0686,xc)));
            x_choque=pto_choque(3,1);
            x_choque_m=pto_choque;
        case 2
            pto_choque=real(double(solve((-b*xc+a*((xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h))==-0.54*(a*xc+b*(xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h).^2+0.0229,xc))); 
            x_choque=pto_choque(2,1);
            x_choque_m=pto_choque;
        case 3
            pto_choque=real(double(solve((-b*xc+a*((xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h))==-0.54*(a*xc+b*(xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h).^2-0.0229,xc))); 
            x_choque=pto_choque(3,1);
            x_choque_m=pto_choque;
        case 4
            pto_choque=real(double(solve((-b*xc+a*((xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h))==-0.54*(a*xc+b*(xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h).^2-0.0686,xc))); 
            x_choque=pto_choque(2,1);
            x_choque_m=pto_choque;
        case 5
            pto_choque=real(double(solve((-b*xc+a*((xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h))==-0.54*(a*xc+b*(xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h).^2-0.1143,xc))); 
            x_choque=pto_choque(3,1);
            x_choque_m=pto_choque;
        case 6
            pto_choque=real(double(solve((-b*xc+a*((xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h))==-0.54*(a*xc+b*(xc-d)*tan(ang_pos)-((9.8*(xc-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h).^2+0.16,xc))); 
            x_choque=pto_choque(2,1);
            x_choque_m=pto_choque;
     end
     y_choque=(x_choque-d)*tan(ang_pos)-((9.8*(x_choque-d)^2)/(2*v^2*(cos(ang_pos)^2)))+h;
end

% yc=-0.54*pto_choque(2,1)^2+0.0686;

% x_choque=pto_choque;

