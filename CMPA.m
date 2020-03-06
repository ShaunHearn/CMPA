function [] = CMPA()
%Diode Parameter Extraction PA session March 6, 2020
%  

%{
    Student:    Shaun Hearn
    Student Id: 100953334

    Course:     ELEC 4700
    Intent:     Diode Parameter Extraction PA
%}

%{ 
 ----------------------------------------------------------------------
                           Constants
------------------------------------------------------------------------
%}
clear
clc
close all
warning off

%{ 
 ----------------------------------------------------------------------
                           Given Values
------------------------------------------------------------------------
%}

Vmin = -1.95;           % minumum voltage (V)
Vmax = 0.7;             % maximum voltage (V)
N = 200;                % number of steps in voltage
variation = 0.2;        % variations in current vector

Is = 0.01e-12;          % saturation current Amps
Ib = 0.1e-12;           % breakdown saturation current Amps
Vb = 1.3;               % breakdown Voltage (V)
Gp = 0.1;               % parasitic parallel capacitance Seimens


%{ 
 ----------------------------------------------------------------------
                           Main code
------------------------------------------------------------------------
%}

parameters = [Is,Gp,Ib,Vb];         % vector of parameters to pass 
V = linspace(Vmin,Vmax,N);          % vector of voltages
I1 = Is*(exp(1.2*V/0.025)-1)+Gp*V-Ib*(exp((-1.2*(V+Vb))/0.025)-1);          % current vector
I2 = I1.*((1-variation/2)+rand(1,length(V))*variation);                     % current vector representing noisey data

numFits = [2,3,1];
for n = 1:3
    figure('units','normalized','outerposition',[0 0 1 1])
    [CurrentplotData,LegendStrings] = generateFits(V,I1,I2,parameters,n,numFits);
    plotFits(V,CurrentplotData,LegendStrings,numFits,n);
end
end


%{ 
 ----------------------------------------------------------------------
                        Function to genreate fits
------------------------------------------------------------------------
%}
function [CurrentplotData,LegendStrings] = generateFits(V,I1,I2,parameters,n,numFits)
LegendStrings = cell(1,numFits(n)+1);
CurrentplotData = cell(numFits(n)+1,2);
for zz = 1:2
    if zz == 1
        IAnalze = I1;
    elseif zz == 2
        IAnalze = I2;
    end
    for xx = 0:numFits(n)
        if xx == 0
            Ifit = IAnalze;
            LegendStrings{1} = 'Data';
        else
            [Ifit,legendStr] = getFitData(V,IAnalze,parameters,n,xx);
            LegendStrings{xx+1} = legendStr;
        end
        CurrentplotData{xx+1,zz} = Ifit;
    end
end
end



%{ 
 ----------------------------------------------------------------------
                        Function to genreate plot fits
------------------------------------------------------------------------
%}
function [] = plotFits(V,CurrentplotData,LegendStrings,numFits,n)
for kk = 1:4
    subplot(2,2,kk)
    for cc = 1:(numFits(n)+1)
        switch kk
            case 1
                plot(V,CurrentplotData{cc,1});
            case 2
                plot(V,CurrentplotData{cc,2});
            case 3
                semilogy(V,abs(CurrentplotData{cc,1}));
            case 4
                semilogy(V,abs(CurrentplotData{cc,2}));
        end
        if cc == 1
           lims = axis; 
        end
        hold on;
    end
    axis(lims);
    xlabel('Voltage (V)')
    ylabel('Current (A)')
    grid on;
    legend(LegendStrings,'Location','southoutside','Orientation','horizontal');
    title(generateTitle(n,kk));
end
end




%{ 
 ----------------------------------------------------------------------
               Function to genreate titles of plots
------------------------------------------------------------------------
%}
function [out] = generateTitle(n,k)
    if n == 1
        out = 'Polynomial Fits';
    elseif n == 2
        out = 'Curve Fits';
    elseif n == 3
        out = 'Neural Network';
    end
    out = strcat(out,': IV Curve ');
    if k == 1  || k == 3
        out = strcat(out,'(No Random Variation)');
    else
        out = strcat(out,'(With Random Variation)');
    end
end


%{ 
 ----------------------------------------------------------------------
               Function to genreate fit options
------------------------------------------------------------------------
%}
function [fo,ft,LegendStrings] = createFitOptions(parameters,mode)
if mode == 1
    upper = [Inf,parameters(2),Inf,parameters(4)];
    lower = [-Inf,parameters(2),-Inf,parameters(4)];
    LegendStrings = 'B and D set';
elseif mode == 2
    upper = [Inf,Inf,Inf,parameters(4)];
    lower = [-Inf,-Inf,-Inf,parameters(4)];
    LegendStrings = 'D set';
else
    upper = [Inf,Inf,Inf,Inf];
    lower = [-Inf,-Inf,-Inf,-Inf];
    LegendStrings = 'No set';
end
ft = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');
fo = fitoptions(ft);
fo = fitoptions(fo,'Lower',lower);
fo = fitoptions(fo,'Upper',upper);
end



%{ 
 ----------------------------------------------------------------------
               Function to access fit data
------------------------------------------------------------------------
%}
function [fitCurrent,LegendStrings] = getFitData(V,I,parameters,questionNum,num)
if questionNum == 1
    p = polyfit(V,I,4*num);
    fitCurrent = polyval(p,V);
    LegendStrings = sprintf('%dth Order Fit',4*num);
elseif questionNum == 2
    [fo,ft,LegendStrings] = createFitOptions(parameters,num);
    ff = fit(V',I',ft,fo);
    fitCurrent = ff(V);
elseif questionNum == 3
    inputs = V;
    targets = I;
    hiddenLayerSize = 10;
    net = fitnet(hiddenLayerSize);
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    net = train(net,inputs,targets);
    outputs = net(V);
    LegendStrings = 'Neural Network Fit';
    fitCurrent = outputs;
end
end