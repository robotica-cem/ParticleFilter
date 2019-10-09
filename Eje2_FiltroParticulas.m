function Eje2_FiltroParticulas
    clc;      close all
    landmarks=[20 80 20 80;20 80 80 20];       mundo=100;

    ruidoMovimiento=[0.2,0.5];       ruidoSenores=5;
    NumParticulas=100;
    % Se crea un robot aleatorio
    myrobot=CreateParticules('Uniforme',[0 mundo;0 mundo; -pi pi],1);
    % Se crean particulas
    P=CreateParticules('Uniforme',[0 mundo;0 mundo; -pi pi],NumParticulas); 
    figure(1); PlotLandMarks(landmarks,mundo,'r');
    PlotParticulas(P,'g'); PlotRobot(myrobot,'k');   pause

    for k=1:100
        d_Medidas=Sensar(myrobot,ruidoSenores,landmarks);
        
        % Corrección
        d_Particulas=Sensar(P,0,landmarks);
        weights=Prob_Medidas(d_Particulas,ruidoSenores,d_Medidas);
        suma=sum(weights);
        P=ReSampling(P,weights/suma,'Sistematico');
        BestEstim=BestRobotPrediccion(P,'promedio');
        
        % Graficar
        figure(1); PlotLandMarks(landmarks,mundo,'r');
        PlotParticulas(P,'g'); PlotRobot(myrobot,'k'); 
        PlotRobot(BestEstim,'r');
        nivel=0.1*NumParticulas*exp(-3^2/2);
        if suma>nivel
            title('LO ENCONTRÉ');
        else
            title('Lo perdí')
        end
        
        pause(0.1)
        
        % Predicción
        Movimiento=[0.1*rand, 5*rand];
        myrobot=Predecir(myrobot,Movimiento,ruidoMovimiento,mundo);
        P=Predecir(P,Movimiento,ruidoMovimiento,mundo);        
    end
end

function best=BestRobotPrediccion(P,metodo)
    switch metodo
        case 'promedio'
            N=size(P,2);
            SumSin=sum(sin(P(3,:)));
            SumCos=sum(cos(P(3,:)));
            best=[mean(P(1,:));mean(P(2,:));...
                  atan2(SumSin/N,SumCos/N)];
        case 'moda'
            best=[mode(P(1,:));mode(P(2,:));mode(P(3,:))];
    end
end


function P=CreateParticules(tipo,Param,numParticulas)
% Esta función crea un grupo de partículas (estado) aleatorias cuya
% dimensión es igual al número de rengloes de la dimensión de la matriz 
% Param.
% 
% numParticulas es el número de partículas que se desean general.
% 
% Si tipo es 'Gaussiana' entonces Param tiene dos columnas, la
%       primera es la media y la segunda es la desviación estandar.
% Si tipo es 'Uniforme' entonces Param tiene dos columnas, la
%       primera contiene el malos mínimo y la segunda el valor máximo.

    numEstados=size(Param,1);
    P=zeros(numEstados,numParticulas);
    switch tipo 
        case'Gaussiana'
            for k=1:numParticulas
                P(:,k)=Param(:,1)+...
                       Param(:,2).*randn(numEstados,1);
            end
        case 'Uniforme'
            for k=1:numParticulas
                P(:,k)=Param(:,1)+(Param(:,2)-...
                        Param(:,1)).*rand(numEstados,1);
            end
        otherwise
            error('No se entiende el tipo de distribución deseada');
    end
end

function P=Predecir(P,u,w,mundo)
% Esta función aplica el modelo de transisción de estado indicado de la
% forma:    P(k+1)=f(P(k), u(k), w(k))
% donde P(k) es el valor del estado para todo el conjunto de partículas
% evaluado en el instante k. Es una matriz de tantas columnas como
% partículas existan.  La acción de control es u.
%
% Ruido w es un vector con la desviación estandar de la incertidumbre del
% movimiento expresado por la ecuación de trasnsición de estado. Este 
% vector es de la misma dimensión que el estado.
%
    numParticulas=size(P,2);
    if u(2)>0
        for k=1:numParticulas
            % Primero gira
            q=P(3,k)+u(1)+w(1)*randn;
            P(3,k)=atan2(sin(q),cos(q));
            % Despues avanza. El mundo es cíclico en ambas direcciones
            a=u(2)+w(2)*randn;
            P(1,k)=EnRango(P(1,k)+a*cos(P(3,k)),mundo);
            P(2,k)=EnRango(P(2,k)+a*sin(P(3,k)),mundo);
        end
    else
        disp('Este robot solo puede avanzar')
    end
end

function medidas=Sensar(P,v,landmarks)
    numLand=size(landmarks,2);
    numParticulas=size(P,2);
    medidas=zeros(numLand,numParticulas);
    for r=1:numParticulas
        for k=1:numLand
            medidas(k,r)=norm(P(1:2,r)-landmarks(:,k))+v*randn;
        end
    end
end

function new_d=EnRango(old_d,ElMax)
    if      old_d>ElMax
        new_d=old_d-ElMax;
    elseif  old_d<0
        new_d=old_d+ElMax;
    else
        new_d=old_d;
    end
end

function PlotLandMarks(landmarks,mundo,color)
    Lx=landmarks(1,:);
    Ly=landmarks(2,:);
    plot(Lx,Ly,[color,'o']);
    axis([-10 mundo+10 -10 mundo+10]);
end

function PlotRobot(robot,color)
    hold on
    plot(robot(1),robot(2),[color,'*'],...
         [robot(1) robot(1)+3*cos(robot(3))],...
         [robot(2) robot(2)+3*sin(robot(3))],[color,'-']);
    hold off
end

function PlotParticulas(P,color)
    hold on
    for k=1:size(P,2)
        plot(P(1,k),P(2,k),[color,'*'],[P(1,k) P(1,k)+3*cos(P(3,k))],...
             [P(2,k) P(2,k)+3*sin(P(3,k))],[color,'-']);
    end
    hold off
end

function weights=Prob_Medidas(d_Particulas,v,d_Medida)
    numParticulas=size(d_Particulas,2);
    weights=zeros(1,numParticulas);
    a=-1/(2*v^2);
    for i=1:numParticulas
        weights(i)=exp(a*norm(d_Particulas(:,i)-d_Medida)^2);
    end
end

function Pnew=ReSampling(Pold,w,Metodo)
    [numEstados,numParticulas]=size(Pold);
    Pnew=zeros(numEstados,numParticulas);
    switch Metodo
        case 'Sistematico2'
            R=cumsum(w,2);
            T=rand(1,numParticulas);
            [~,idx]=histc(T,R);
            Pnew=Pold(:,idx+1);
        case 'Sistematico'
            edges=[0 cumsum(w)]*numParticulas;
            [~, idx] = histc(rand:numParticulas,edges);
            Pnew = Pold(:,idx);
        case 'Multinominal'
            idx = randsample(1:numParticulas,numParticulas,true,w);
            Pnew = Pold(:,idx);
        case 'Acumulado'
            acumul=cumsum(w);
            for i=1:numParticulas
                Pnew(:,i)=Pold(:,find(rand <= acumul,1));
            end
        case 'Ruleta'
            index=1+round((numParticulas-1)*rand);
            beta=0;
            mw=max(w);
            for i=1:numParticulas
                beta=beta+2*mw*rand;
                while beta > w(index)
                    beta=beta-w(index);
                    index=EnRango(index+1,numParticulas);
                end
                Pnew(:,i)=Pold(:,index);
            end
        otherwise
          error('Estrategia de re-muestreo no implementada')
    end
end