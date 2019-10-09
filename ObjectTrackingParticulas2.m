function ObjectTrackingParticulas2
clc
close all;


base_dir = '.\hexbug_frames_compressed\';
cd(base_dir);
f_list=dir('*png'); % obtiene la lista de las im�genes png en la carpeta.
S_frame=10; %Obtendremos una imagen promedio de estos frames
[ren,col,cap]=size(imread(f_list(1).name));
imagen=uint16(zeros(ren,col,cap));
for k=1:S_frame
    laimagen=uint16(imread(f_list(k).name));
    imagen = imadd(imagen,laimagen,'uint16');
end
backImagen=uint8(imagen/S_frame);

% Inicializaci�n de las part�culas
numParticles=400;
desvW=4; desvV=6; % ruido movimiento y sensores
P=NewRobot('Gaussiana',[1 desvW;1 desvW;0 desvW;0 desvW],numParticles);
% valor incial del estado estimado

for k=S_frame+1:length(f_list)
    % Esta secci�n lee la posici�n de insecto y produce su posici�n [x,y]
    imagen=imread(f_list(k).name);
    diffImagen=imabsdiff(backImagen,imagen);
    BW = im2bw(diffImagen,0.15);
    s=regionprops(BW,'centroid','area');
    [~,donde]=max(cat(1,s.Area));
    
    if ~isempty(donde)
        x=s(donde).Centroid(1);
        y=s(donde).Centroid(2);
    else
        x=NaN;
        y=NaN;
    end
    if (k>240)&&(k<280)
        Y=[NaN;NaN];
        imshow(backImagen)
    else
        Y=[x; y];
        imshow(imagen)
    end
    if k==S_frame+1, pause; end
    P=filtroParticulas(P,Y,[col;ren],[desvV;desvW]);  
end
cd ..
end

function Pk=filtroParticulas(Pk,yk,Resolucion,ruido)
    metodoReMuestreo='Sistematico'; % 'Multinominal'; 'Acumulado'; 'Sistematico'; 'Ruleta';
    metodoMejor='Media';  %'Media';  'Moda';  'Centro';
    desvRuidoSensado=ruido(1);
    desvRuidoMov=ruido(2);

    %Correcci�n
    if ~isnan(yk)
        yEstimada=Sensar(Pk,[0;0]);
        weights=Prob_Medidas(yEstimada,yk,desvRuidoSensado);
        Pk=ReSampling(Pk,weights,metodoReMuestreo);
    end
    %Predicci�n
    Pk=Mover(Pk,desvRuidoMov,Resolucion);
    xkBest=BestPrediccion(Pk,metodoMejor);
    showParticules(Pk,xkBest,yk)
end

function showParticules(P,PBest,Y)
    hold on; 
    plot(P(1,:),P(2,:),'g.',PBest(1),PBest(2),'ro',Y(1),Y(2),'k*');
    hold off
    drawnow
end

function P=Mover(P,w,Resolucion)
% Esta funci�n aplica el modelo de transisci�n de estado indicado de la
% forma:    P(k+1)=f(P(k), w(k))
% donde P(k) es el valor del estado para todo el conjunto de part�culas
% evaluado en el instante k. Es una matriz de tantas columnas como
% part�culas existan. 
%
% Ruido w es un vector con la desviaci�n estandar de la incertidumbre del
% movimiento expresado por la ecuaci�n de trasnsici�n de estado. Este 
% vector es de la misma dimensi�n que el estado.
%
    numParticulas=size(P,2);
    T=1;  A=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];  
    for i=1:numParticulas
        P(:,i)=A*P(:,i);
        P(1,i)=1+EnRango(P(1,i)+w*randn,Resolucion(1)-1);
        P(2,i)=1+EnRango(P(2,i)+w*randn,Resolucion(2)-1);
    end
end

function y=Sensar(P,v)
% Esta funci�n produce la se�anl de salida que entegan los sensores
% correspondiente a cada part�cula, cuya forma es:
%  y(k)=g(P(k), v(k))
% donde: P(k) es el conjunto de part�culas evaluado en el instante k. Es 
% una matriz de tantas columnas como part�culas existan.
% Ruido v es un vector con la desviaci�n estandar de la incertidumbre en
% los sensores. Este vector es de la misma dimensi�n que el n�mero de
% sensores.
%
    numParticulas=size(P,2);
    numSensores=size(v,1);
    y=zeros(numSensores,numParticulas);
    for i=1:numParticulas
        y(:,i)= P(1:2,i);
    end
end

function P=NewRobot(tipo,ValoresIniciales,numParticulas)
% Esta funci�n crea un grupo de part�culas (estado) aleatorias cuya
% dimensi�n es igual al n�mero de rengloes de la dimensi�n de la matriz 
% ValoresInciales.
% 
% numParticulas es el n�mero de part�culas que se desean general.
% 
% Si tipo es 'Gaussiana' entonces ValoresIniciales tiene dos columnas, la
%       primera es la media y la segunda es la desviaci�n estandar.
% Si tipo es 'Uniforme' entonces Valores Inciales tiene dos columnas, la
%       primera contiene el malos m�nimo y la segunda el valor m�ximo.

    numEstados=size(ValoresIniciales,1);
    P=zeros(numEstados,numParticulas);
    if strcmp(tipo,'Gaussiana')
        for k=1:numParticulas
            P(:,k)=ValoresIniciales(:,1)+...
                   ValoresIniciales(:,2).*randn(numEstados,1);
        end
    elseif strcmp(tipo,'Uniforme')
        for k=1:numParticulas
            P(:,k)=ValoresIniciales(:,1)+(ValoresIniciales(:,2)-...
                    ValoresIniciales(:,1)).*rand(numEstados,1);
        end
    else
        error('No se entiende el tipo de distribuci�n deseada');
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
          error('Metodo de selecci�n del mejor no implementado')
    end            
end

%
function Pnew=ReSampling(Pold,w,Metodo)
    [numEstados,numParticulas]=size(Pold);
    Pnew=zeros(numEstados,numParticulas);
    switch Metodo
        case 'Sistematico'
            edges = min([0 cumsum(w)],1);
            edges(end) = 1;
            u1 = rand/numParticulas;
            [~, idx] = histc(u1:1/numParticulas:1,edges);
            Pnew = Pold(:,idx);
        case 'Multinominal'
            idx = randsample(1:numParticulas,numParticulas,true,w);
            Pnew = Pold(:,idx);
        case 'Acumulado'
            for i=1:numParticulas
                Pnew(:,i)=Pold(:,find(rand <= cumsum(w),1));
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

function weights=Prob_Medidas(yEstimados,yReales,ruido)
    [numSensores,numParticulas]=size(yEstimados);
    weights=ones(1,numParticulas);
    for i=1:numParticulas
        for k=1:numSensores
            weights(i)=weights(i)*...
                       normpdf(yEstimados(k,i)-yReales(k),0,ruido);
        end
    end
    weights=weights/sum(weights);
end