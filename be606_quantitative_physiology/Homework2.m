%Liam Jackson 
%U40546227
%HW 2, BE606

clear all 
close all 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1c plot

syms t tau1 tau2 pmax taurise;

B = (((tau2/tau1) ^ (taurise/tau1)) - ((tau2/tau1) ^ (taurise/tau2)))^-1;
psb = pmax * B * (exp(-t/tau1)-exp(-t/tau2)); 
psbPlot = matlabFunction(psb);
tmaxb = solve(0==diff(psb)); 
tmaxbF = matlabFunction(tmaxb); %already calculated the tmaxes on paper but want to verify

syms taus;

psa = (pmax * t / taus)*exp(-t/taus);
psaPlot = matlabFunction(psa);
tmaxa = solve(0==diff(psa));
tmaxaF = matlabFunction(tmaxa); %already calculated the tmaxes on paper but want to verify

clear t tau1 tau2 pmax taurise taus;

t = linspace(0, 10, 100);
tau1 = 5.6;
syms tau2dumb; %only for getting a precise result for tau2. On paper I rounded, here I'm going for precision
tau2 = double(solve(0.3 == (tau1*tau2dumb)/(tau1-tau2dumb), tau2dumb));
taurise = (tau1*tau2)/(tau1-tau2);
taus = taurise*log(tau1/tau2);
B = (((tau2/tau1) ^ (taurise/tau1)) - ((tau2/tau1) ^ (taurise/tau2)))^-1;
pmax = 1;

disp(['tmaxA = ', num2str(tmaxaF(taus))])         %verifying if the tmaxes 
disp(['tmaxB = ', num2str(tmaxbF(tau1,tau2))])    %are the same

figure()

plot(t,psbPlot(pmax,t,tau1,tau2,taurise))
hold on 
plot(t, psaPlot(pmax, t,taus))
xline(tmaxaF(taus),'--');
xlim([0 10])
title({'GABA(A) Receptor Open Probability','and Synaptic Conductance vs. time'})
legend('Open Probability (Receptor)','Synaptic Conductance Probability', 't_{max} = 0.894 sec')
ylabel('Probability (arb)')
xlabel('Time (ms)')
grid on 
hold off

%With both curves normalized to Pmax, PsB reaches Pmax and has a much 
%longer decay. PsA does not reach Pmax, and has a sharper decay.
%This could be due to concentration gradients working against the
%conductance even though the Receptors are fully open

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%2a.

%Parameters
Cm = 10E-3; 
Ie = 0;
gl = 19E-3; 
El = -67E-3; 
gna = 74E-3;
Ena = 60E-3; 
Vhalf = 1.5E-3;
k = 16E-3;
tspan = [0 20];
IC1 = [.05; -10]; 

%Symbolic Variables
syms Em(t) 

%Initializing arrays
minf = [];
Il = [];
Inap = [];
Im = [];

%Solve ODE
[tspan, Em] = ode45(@odefxn, tspan, IC1);
Em = Em(:,1);

for ind = 1:length(Em(:,1))
    minf(ind) = 1 / (1+exp((Vhalf-Em(ind)) / k));
    Il(ind) = gl*(Em(ind) - El);
    Inap(ind) = gna * minf(ind) * (Em(ind)-Ena);
    Im(ind) = Il(ind) + Inap(ind);
end

%Plotting
figure()
plot(Em, Il)
hold on 
plot(Em, Inap)
plot(Em, Im)
xline(-.04132, '--');
xline(.00289, '--');
title('Current Trends vs. Membrane Voltage')
legend('Il','Inap','Im', 'Region of Pos Feedback')
ylabel('Current (A)')
xlabel('Voltage (V)')
grid on
hold off

%The slop of the 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2b.

figure()
IC2 = [-.017; -.000418];    %these ICs give a good representation of the 
                            %AP Upstroke and were found in the neg slope 
                            %region of Im in the preceding figure
%Initializing arrays
minf2 = [];
Il2 = [];
Inap2 = [];
Im2 = [];

%Solve ODE
[tspan, Em2] = ode45(@odefxn, tspan, IC2);
Em2 = Em2(:,1);

for ind = 1:length(Em2(:,1))
    minf2(ind) = 1 / (1+exp((Vhalf-Em2(ind)) / k));
    Il2(ind) = gl*(Em2(ind) - El);
    Inap2(ind) = gna * minf2(ind) * (Em2(ind)-Ena);
    Im2(ind) = Il2(ind) + Inap2(ind);
end

plot(tspan, Em2)
hold on 
plot(tspan, Im2)
title('Membrane Voltage and Current vs. Time')
legend('Em','Im')
xlim([-0.5 20])
xlabel('Time (ms)')
ylabel('Intensity')
grid on
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%3a.
SimLen = 10;     % Length of simulation (seconds)
dt = .001;      % 1ms resolution
r0 = 100;        % spike frequency (Hz)
nTrials = 10;   % number of trials 

%Generate Spike Trains 
[SpikeMatrix, TotalTimeVector] = SpikeGen(r0,SimLen,nTrials);

%Find Spike Occurrences
for i = 1:nTrials
    for j = 1:length(SpikeMatrix(1,:))
        [Trial, TimeofSpike] = find(SpikeMatrix == 1);
    end
end

[Trial,sortIndex] = sort(Trial,'ascend');
TimeofSpike = TimeofSpike(sortIndex);

%Spike times for each trial 
SpikeTimesPerTrial = [Trial, TimeofSpike];

%Calc Interspike Intervals
IntervalsRaw = diff(SpikeTimesPerTrial(:,2));
InterspikeIntervals = IntervalsRaw(IntervalsRaw > 0);

%Plot Histogram
figure()
ISI = histogram(InterspikeIntervals);
title({'ISI Distribution of Neuron with','Constant Firing Rate','(No Refractoriness)'})
xlabel('Interval Length (ms)')
ylabel('Number of Occurrences')

%Raster Plot
figure() 
Raster(SpikeMatrix, TotalTimeVector);
title({'Raster Plot of Neuron with','Constant Firing Rate','(No Refractoriness)'})
xlabel('Time (sec)')
ylabel('Trial')
CoeffVar = std(InterspikeIntervals)/mean(InterspikeIntervals)

%3b.

%These are the same parameters as before:
% SimLen = 10;         % Length of simulation (seconds)
% dt = .001;          % 1ms resolution
% r0 = 100;           % initial fire rate
tauref = .01;         % Refractory Recovery Rate (seconds)
% nTrials = 10;       % number of trials 

%Generate Spike Trains 
[SpikeMatrixRef, TotalTimeVector] = SpikeGenRef(r0, tauref, SimLen, dt, nTrials);

%Find Spike Occurrences
for i = 1:nTrials
    for j = 1:length(SpikeMatrixRef(1,:))
        [Trial, TimeofSpike] = find(SpikeMatrixRef == 1);
    end
end

[Trial,sortIndex] = sort(Trial,'ascend');
TimeofSpike = TimeofSpike(sortIndex);

%Spike times for each trial 
SpikeTimesPerTrialRef = [Trial, TimeofSpike];

%Calculate Interspike Intervals
IntervalsRawRef = diff(SpikeTimesPerTrialRef(:,2));
InterspikeIntervalsRef = IntervalsRawRef(IntervalsRawRef > 0);

%Plot Histogram
figure()
ISIRef = histogram(InterspikeIntervalsRef, 50);
title({'ISI of Neuron with','Time Dependent Firing Rate','(Refractory Period)'})
xlabel('Interval Length (ms)')
ylabel('Number of Occurrences')

%Raster Plot
figure()
RasterRef(SpikeMatrixRef, TotalTimeVector);
title({'Raster Plot of Neuron with','Time Dependent Firing Rate','(Refractory Period)'})
xlabel('Time (sec)')
ylabel('Trial')
CoeffVarRef = std(InterspikeIntervalsRef)/mean(InterspikeIntervalsRef)

%Panel Plot just because it's interesting to compare
figure()
subplot(2,2,1)
histogram(InterspikeIntervals);
title({'ISI Distribution of Neuron with','Constant Firing Rate','(No Refractoriness)'})
xlabel('Interval Length (ms)')
ylabel('Number of Occurrences')

subplot(2,2,3)
Raster(SpikeMatrix, TotalTimeVector);
title({'Raster Plot of Neuron with','Constant Firing Rate','(No Refractoriness)'})
xlabel('Time (sec)')
ylabel('Trial')

subplot(2,2,2)
histogram(InterspikeIntervalsRef, 50);
title({'ISI Distribution of Neuron with','Time Dependent Firing Rate','(Refractory Period)'})
xlabel('Interval Length (ms)')
ylabel('Number of Occurrences')

%The refractory period changes the distribution of the histogram, causing a
%skew to the right. Also, the addition of a refractory period greatly 
%increases the average length of interval between spikes.

subplot(2,2,4)
RasterRef(SpikeMatrixRef, TotalTimeVector);
title({'Raster Plot of Neuron with','Time Dependent Firing Rate','(Refractory Period)'})
xlabel('Time (sec)')
ylabel('Trial')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%4
RawData = importdata('C1P8data.xlsx');

StimSpikeArray = [RawData(:,1) RawData(:,2)]; %[stim spikes]

WindowSize = 150; %150 time steps ie 300ms windows
SampleRate = 500; %Hz

%Find Spike Occurrences
AllSpikeIndices = find(StimSpikeArray(:,2) == 1); %col vect of spike indices from rawdata 

%Initialize Matrix for Storing Windows
WindowSum = zeros(WindowSize,1);

%Sum stimulus windows prior to each spike
for RowOfStimSpikeArray = WindowSize:length(StimSpikeArray(:,2)) %row 150 thru 60k
    if StimSpikeArray(RowOfStimSpikeArray,2) == 1
        WindowSum = WindowSum + StimSpikeArray((RowOfStimSpikeArray-WindowSize:RowOfStimSpikeArray-1),1);
    end
end

%Average Response-Inducing Stimulus = sum of all windows / number of spikes
AvgStim = WindowSum ./ length(AllSpikeIndices);

WindowLength = (1000 * WindowSize / SampleRate); % (1000ms/1sec)*(150 samples/500 samples/sec) 
WindowLengthMap = 1:2:WindowLength; %time vector for plotting

figure()
plot(WindowLengthMap, AvgStim)
title({'Spike-Triggered Average','Response to Random Stimulus'})
xlabel('Stimulus Window Duration (milliseconds)')
ylabel('Stimulus Intensity')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Functions

function [SpikeMatrix, TimeVec] = SpikeGen(r, SimLen, nTrials)
    dt = 1/1000; % 1 ms resolution
    nBins = floor(SimLen/dt);
    SpikeMatrix = rand(nTrials, nBins) < r*dt;
    TimeVec = 0:dt:SimLen-dt;
end

function [SpikeMatrixRef, TimeVec] = SpikeGenRef(r0, tauref, SimLen, dt, nTrials)
    nBins = floor(SimLen/dt);
    TimeVec = 0:dt:SimLen-dt;
    t = TimeVec;
    r = zeros(1,length(t)); 
    for t = 2:(length(t)-1)
        r(t+1) =  r0 + (r(t)-r0)*exp(-dt/tauref);
        for y = 1:nTrials
            for x = 1:length(TimeVec)
                SpikeMatrixRef(y,x) = rand; 
                if SpikeMatrixRef(y,x) < r(t) * dt
                    SpikeMatrixRef(y,x) = 1;
                    r(t+1) = 0;
%                     t = t+1;
                else 
                    SpikeMatrixRef(y,x) = 0;
%                     t = t+1;
                end
            end
        end
    end
end

function [] = Raster(SpikeMatrix, TimeLen)
hold all;
for trialCount = 1:size(SpikeMatrix,1)
    spikePos = TimeLen(SpikeMatrix(trialCount,:) == 1);
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0.5 size(SpikeMatrix, 1)+.5]);
end

function [] = RasterRef(SpikeMatrixRef, TimeLen)
hold all;
for trialCount = 1:size(SpikeMatrixRef,1)
    spikePos = TimeLen(SpikeMatrixRef(trialCount,:) == 1);
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0.5 size(SpikeMatrixRef, 1)+.5]);
end

function dEmdt = odefxn(tspan,Em)
%Parameters
Cm = 10E-3; 
Ie = 0;
gl = 19E-3; 
El = -67E-3; 
gna = 74E-3;
Ena = 60E-3; 
Vhalf = 1.5E-3;
k = 16E-3;

dEmdt = (-1/Cm)*(gl * (Em - El) + ...
    (gna * (1/(1+exp((Vhalf-Em)/k))) * (Em - Ena)));
end
