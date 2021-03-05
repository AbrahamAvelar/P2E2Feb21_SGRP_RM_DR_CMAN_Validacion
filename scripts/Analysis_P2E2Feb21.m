for platos=[1:4 6:7]
    y=tic;
    %close all
    directorio='D:\Dropbox\Posdoc\Experimentos\P2E2Feb21_SGRP_RM_DR_CMAN_Validacion\E2Feb21_SGRP_RM_DR_AND_CMANValidation\Experiment_068';
    PL = LoadCitPlate(directorio, platos);
    platename=strcat('PL',num2str(platos))
    cd 'D:\Dropbox\Posdoc\Experimentos\P2E2Feb21_SGRP_RM_DR_CMAN_Validacion'
    save (platename,'-v7.3')
    toc(y)
end

%% Consolidar datos en AllData
con=0;
for plato={'PL1','PL2','PL3','PL4','PL6','PL7'}
	load((cell2mat(plato)))
    con=con+1;
    AllData(con).PL=PL;
end
save ('AllData','AllData','-v7.3')
% Quitar valores negativos en x & y
y=3; %SYBR Green
x=4; %PI
AllDataNoLog = QuitaLogNeg(AllData, x, y);
save ('AllDataNoLog','AllDataNoLog','-v7.3')
%%  Hacer Gates para todas las réplicas
y=7; %SYBR Green
x=8; %PI
namex=AllDataNoLog(1).PL(1).Info.par(x).name;
namey=AllDataNoLog(1).PL(1).Info.par(y).name;

[GateArrayV, GateArrayM]=TwoGatesSubset(AllDataNoLog(1).PL(3), wells, samplesize, x, y, namex, namey);


%% Calcular los porcentajes con los gates de la sección pasada
y=7; %SYBR Green
x=8; %PI
namex=AllDataNoLog(1).PL(1).Info.par(x).name;
namey=AllDataNoLog(1).PL(1).Info.par(y).name;
clear AllData
for plato=1:length(AllDataNoLog)
    for pl=1:length(AllDataNoLog(plato).PL)
        for w=1:size(AllDataNoLog(plato).PL(pl).WELL,2)
        columna=1;
        [GatedDataV, GatedIndexesV] = ApplyGate(log10(AllDataNoLog(plato).PL(pl).WELL(w).dat), GateArrayV, x, y);
        [GatedDataM, GatedIndexesM] = ApplyGate(log10(AllDataNoLog(plato).PL(pl).WELL(w).dat), GateArrayM, x, y);
        
        V=length(GatedIndexesV);
        M=length(GatedIndexesM);
        AllData(plato).pctV(pl,w) = V/(V+M);
        AllData(plato).pctM(pl,w) = M/(V+M);
        
        AllData(plato).events(pl,w) = size(log10(AllDataNoLog(plato).PL(pl).WELL(w).dat),1);
        AllData(plato).pctVtot(pl, w) = V/size(log10(AllDataNoLog(plato).PL(pl).WELL(w).dat),1);
        AllData(plato).pctMtot(pl, w) = M/size(log10(AllDataNoLog(plato).PL(pl).WELL(w).dat),1);
        end
    end
end

GatedPctgs = AllData;
save GatedPctgs GatedPctgs

%% sacar los vectores de tiempo
for plate = 1:length(AllDataNoLog)
    for pl = 1: length(AllDataNoLog(plate).PL)
        temp=strsplit(AllDataNoLog(plate).PL(pl).Info.PlateName,'_');
        GatedPctgs(plate).t(pl) = datenum(cell2mat(temp(2)), 'YYYYMMDD');
    end
    GatedPctgs(plate).t=GatedPctgs(plate).t- GatedPctgs(plate).t(1);
end

%% curvas de muerte
for pl=1:6
    figure(pl)
    clf
    for w=1:96
        subplot(8,12,w)
        plot(GatedPctgs(pl).t, GatedPctgs(pl).pctV(:,w))
        ylim([.5 1.])
        xlim([0 2])
        titulo=strsplit(AllDataNoLog(pl).PL(1).WELL(w).info.filename, '_');
        title(titulo(3))
    end
end

%% Número de eventos
for pl=6
    figure(pl+10)
    clf
    %for dia=1:size(AllData(pl).events,1)
        for w=1:96
            subplot(8,12,w)
            plot(AllData(pl).events(:,w))
            ylim([0 20000])
            titulo=strsplit(AllDataNoLog(pl).PL(1).WELL(w).info.filename, '_');
            set(gca, 'yticklabel',[])
            title(titulo(3))
        end
    %end
end

%% ver gates de algo específico
AllDataNoLog(1).GateArrays((1)).vivas=GateArrayV;
AllDataNoLog(1).GateArrays((1)).muertas=GateArrayM;
% ver gates de algo específico
platos=1;
pls = 1:3;
%ws=[1 4 7 10];
%ws=[3 6 9 12];
ws=87:3:96;
samplesize=1000;
PlotGatedScatter(AllDataNoLog, platos, pls, ws, x, y, samplesize);


