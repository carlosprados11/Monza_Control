%%

load('dificultad1.mat');
load('circulos.mat');
giro_simu=giro_m;
ma=cos(giro_simu);
mb=sin(giro_simu);

anim_giro_pistas;


%**************************************************************************
%***************GR�FICAS DE LAS CURVAS DE LA SIMULACI�N********************
%**************************************************************************

fig=figure(1);
for j=1:length(xs)
    plot(mxp(j),myp(j)+0.01,'r*','lineWidth',15)
    
    %*********RECTAS*********
    plot(mxl1(j,:),myl1(j,:),'b','lineWidth',2);
    plot(mxl2(j,:),myl2(j,:),'b','lineWidth',2);
    plot(mxl3(j,:),myl3(j,:),'b','lineWidth',2);
    plot(mxl4(j,:),myl4(j,:),'b','lineWidth',2); 
    %*******PARABOLAS***********
    plot(mx(j,:),my(j,:),'b','lineWidth',2);
    plot(mx1(j,:),my1(j,:),'b','lineWidth',2);
    plot(mx2(j,:),my2(j,:),'b','lineWidth',2);
    plot(mx3(j,:),my3(j,:),'b','lineWidth',2);
    plot(mx4(j,:),my4(j,:),'b','lineWidth',2);
    plot(mx5(j,:),my5(j,:),'b','lineWidth',2);
    plot(mx6(j,:),my6(j,:),'b','lineWidth',2);
    plot(mx7(j,:),my7(j,:),'b','lineWidth',2);
      
    hold off;
    MM = getframe(fig);
    figure(1)
    plot(r1,r2,'r');
    hold on;
    plot(r3,r4,'r');
    hold on;
    axis([-0.25 0.25 -0.25 0.25])
    title('SIMULADOR GR�FICO PLATAFORMA PARA PRUEBA Y ENSAYO DE CONTROLADORES')
end

