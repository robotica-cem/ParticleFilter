function ObjectTrackingParticulas
close all;

vr = VideoReader('.\PF_Video_EN\Person.wmv');
col=vr.Width;  ren=vr.Height;
Nfrm_movie = floor(vr.Duration * vr.FrameRate);

% Inicialización de las partículas
desvWpos=25; desvWvel=5; desvV=50; % ruido movimiento y sensores
numParticulas=20000;  ColorDeseado=[255;0;0];
P=CreateParticules('Uniforme',[1 col;1 ren;0 0;0 0],numParticulas);

tic
for k = 1:Nfrm_movie
    
    Frame = readFrame(vr,'default'); 
    % Filtro de Particulas
    [weights,suma]=Prob_Medidas(P,Frame,ColorDeseado,desvV,[col;ren]);
    P=ReSampling(P,weights/suma,'Sistematico');
    % Dibujar Particulas
    showParticles(P,Frame,suma);
    P=Predecir(P,[desvWpos;desvWvel]);  
    if k==1, pause; end
end
toc
end

function showParticles(P,Frame,suma)
    image(Frame);
    metodoMejor='Media';
    bestP=BestPrediccion(P,metodoMejor);
    hold on; 
    plot(P(1,:),P(2,:),'g.',bestP(1),bestP(2),'k+'); 
    hold off;
    numParticulas=size(P,2);
    nivel=0.1*numParticulas*exp(-3^2/2);
    if suma>nivel
        title('Objeto encontrado')
    else
        title('El objeto buscado no está');
    end
    drawnow
end

function P=Predecir(P,w)
% Esta función aplica el modelo de transición de estado indicado de la
% forma:    P(k+1)=f(P(k), w(k))
% donde P(k) es el valor del estado para todo el conjunto de partículas
% evaluado en el instante k. Es una matriz de tantas columnas como
% partículas existan. 
%
% Ruido w es un vector con la desviación estandar de la incertidumbre del
% movimiento expresado por la ecuación de trasnsición de estado. Este 
% vector es de la misma dimensión que el estado.
%
    numParticulas=size(P,2);
    A=[1 0 1 0;0 1 0 1;0 0 1 0;0 0 0 1];
    P=A*P;
    P(1:2,:)=ceil(P(1:2,:)+w(1)*randn(2,numParticulas));
    P(3:4,:)=P(3:4,:)+w(2)*randn(2,numParticulas);
end

function P=CreateParticules(tipo,ValoresIniciales,numParticulas)
% Esta función crea un grupo de partículas (estado) aleatorias cuya
% dimensión es igual al número de rengloes de la dimensión de la matriz 
% ValoresInciales.
% 
% numParticulas es el número de partículas que se desean general.
% 
% Si tipo es 'Gaussiana' entonces ValoresIniciales tiene dos columnas, la
%       primera es la media y la segunda es la desviación estandar.
% Si tipo es 'Uniforme' entonces Valores Inciales tiene dos columnas, la
%       primera contiene el malos mínimo y la segunda el valor máximo.

    numEstados=size(ValoresIniciales,1);
    P=zeros(numEstados,numParticulas);
    switch tipo 
        case'Gaussiana'
            for k=1:numParticulas
                P(:,k)=ValoresIniciales(:,1)+...
                       ValoresIniciales(:,2).*randn(numEstados,1);
            end
        case 'Uniforme'
            for k=1:numParticulas
                P(:,k)=round(ValoresIniciales(:,1)+(ValoresIniciales(:,2)-...
                        ValoresIniciales(:,1)).*rand(numEstados,1));
            end
        otherwise
            error('No se entiende el tipo de distribución deseada');
    end
end

function bestParticula=BestPrediccion(P,Metodo)
    switch Metodo
        case 'Media'
            bestParticula=mean(P(1:2,:),2);
        case 'Moda'
            bestParticula=mode(P(1:2,:),2);
        case 'Centro'
            bestParticula=(min(P(1:2,:),[],2)+max(P(1:2,:),[],2))/2;
       otherwise
          error('Metodo de selección del mejor no implementado')
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

function new_d=EnRango(old_d,ElMax)
    if      old_d>ElMax
        new_d=old_d-ElMax;
    elseif  old_d<0
        new_d=old_d+ElMax;
    else
        new_d=old_d;
    end
end

function [weights,suma]=Prob_Medidas(P,Frame,ColorDeseado,ruido,Res)
    numParticulas=size(P,2);
    weights=ones(1,numParticulas);
    %a=1/(sqrt(2*pi)*ruido);
    b=-0.5/(ruido^2);
    for i=1:numParticulas
        r=P(1,i);
        c=P(2,i);
        if (r>=1 && r<=Res(1)) && (c>=1 && c<=Res(2))
            y=double([Frame(c,r,1);Frame(c,r,2);Frame(c,r,3)]); %Sensar
            d2=(y(1)-ColorDeseado(1))^2+...
               (y(2)-ColorDeseado(2))^2+...
               (y(3)-ColorDeseado(3))^2;
            weights(i)=exp(b*d2);
        else
            weights(i)=0;
        end
    end
    suma=sum(weights);
end