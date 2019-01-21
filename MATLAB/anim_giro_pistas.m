%%**************************ANIMACIÓN DE LAS PISTAS GIRANDO****************

%***********RECTAS******************
mxl1=zeros(length(giro_simu),length(xl1));
myl1=zeros(length(giro_simu),length(yl1));
for i=1:length(giro_simu)
    for j=1:length(xl1)
        mxl1(i,j)=xl1(j)*ma(i)-yl1(j)*mb(i);
        myl1(i,j)=xl1(j)*mb(i)+yl1(j)*ma(i);
    end
end

mxl2=zeros(length(giro_simu),length(xl2));
myl2=zeros(length(giro_simu),length(yl2));
for i=1:length(giro_simu)
    for j=1:length(xl2)
        mxl2(i,j)=xl2(j)*ma(i)-yl2(j)*mb(i);
        myl2(i,j)=xl2(j)*mb(i)+yl2(j)*ma(i);
    end
end

mxl3=zeros(length(giro_simu),length(xl3));
myl3=zeros(length(giro_simu),length(yl3));
for i=1:length(giro_simu)
    for j=1:length(xl3)
        mxl3(i,j)=xl3(j)*ma(i)-yl3(j)*mb(i);
        myl3(i,j)=xl3(j)*mb(i)+yl3(j)*ma(i);
    end
end

mxl4=zeros(length(giro_simu),length(xl4));
myl4=zeros(length(giro_simu),length(yl4));
for i=1:length(giro_simu)
    for j=1:length(xl4)
        mxl4(i,j)=xl4(j)*ma(i)-yl4(j)*mb(i);
        myl4(i,j)=xl4(j)*mb(i)+yl4(j)*ma(i);
    end
end


%********PARABOLAS****************
mx=zeros(length(giro_simu),length(xp));
my=zeros(length(giro_simu),length(yp));
for i=1:length(giro_simu)
    for j=1:length(xp)
        mx(i,j)=xp(j)*ma(i)-yp(j)*mb(i);
        my(i,j)=xp(j)*mb(i)+yp(j)*ma(i);
    end
end


mx1=zeros(length(giro_simu),length(xp1));
my1=zeros(length(giro_simu),length(yp1));
for i=1:length(giro_simu)
    for j=1:length(xp1)
        mx1(i,j)=xp1(j)*ma(i)-yp1(j)*mb(i);
        my1(i,j)=xp1(j)*mb(i)+yp1(j)*ma(i);
    end
end

mx2=zeros(length(giro_simu),length(xp2));
my2=zeros(length(giro_simu),length(yp2));
for i=1:length(giro_simu)
    for j=1:length(xp2)
        mx2(i,j)=xp2(j)*ma(i)-yp2(j)*mb(i);
        my2(i,j)=xp2(j)*mb(i)+yp2(j)*ma(i);
    end
end

mx3=zeros(length(giro_simu),length(xp3));
my3=zeros(length(giro_simu),length(yp3));
for i=1:length(giro_simu)
    for j=1:length(xp3)
        mx3(i,j)=xp3(j)*ma(i)-yp3(j)*mb(i);
        my3(i,j)=xp3(j)*mb(i)+yp3(j)*ma(i);
    end
end

mx4=zeros(length(giro_simu),length(xp4));
my4=zeros(length(giro_simu),length(yp4));
for i=1:length(giro_simu)
    for j=1:length(xp4)
        mx4(i,j)=xp4(j)*ma(i)-yp4(j)*mb(i);
        my4(i,j)=xp4(j)*mb(i)+yp4(j)*ma(i);
    end
end

mx5=zeros(length(giro_simu),length(xp5));
my5=zeros(length(giro_simu),length(yp5));
for i=1:length(giro_simu)
    for j=1:length(xp5)
        mx5(i,j)=xp5(j)*ma(i)-yp5(j)*mb(i);
        my5(i,j)=xp5(j)*mb(i)+yp5(j)*ma(i);
    end
end

mx6=zeros(length(giro_simu),length(xp6));
my6=zeros(length(giro_simu),length(yp6));
for i=1:length(giro_simu)
    for j=1:length(xp6)
        mx6(i,j)=xp6(j)*ma(i)-yp6(j)*mb(i);
        my6(i,j)=xp6(j)*mb(i)+yp6(j)*ma(i);
    end
end

mx7=zeros(length(giro_simu),length(xp7));
my7=zeros(length(giro_simu),length(yp7));
for i=1:length(giro_simu)
    for j=1:length(xp7)
        mx7(i,j)=xp7(j)*ma(i)-yp7(j)*mb(i);
        my7(i,j)=xp7(j)*mb(i)+yp7(j)*ma(i);
    end
end

%***************ROTACION DE LA MONEDA**************************

mxp=zeros(1,length(xs));
for i=1:length(xs)
%     if flag_mov(i)==0
%       mxp(i)=xs(i)*ma(i)-ys(i)*mb(i);
%     else
         mxp(i)=xs(i);
%     end
end
myp=zeros(1,length(ys));
for i=1:length(ys)
%     if flag_mov(i)==0
%       myp(i)=xs(i)*mb(i)+ys(i)*ma(i);
%     else
        myp(i)=ys(i);
%     end
end