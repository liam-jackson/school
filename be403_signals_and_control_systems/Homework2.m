%Question 3c
load pulse_wave(1).mat
p = pulse_wave;
t = linspace(0,(length(p)-1)/300,length(p));
H = zeros(length(p));
window = 3; 
H([1:window],1) = 1;

for i = 2:length(p)
    Hi = H(:,i-1);
    Hi = circshift(Hi,1);
    Hi(1) = 0;
    H(:,i) = Hi;
end

pfiltered = (1/window)*(H*p);

figure
plot(t,p)
hold on 
plot(t,pfiltered)
title('Pulsewave Raw Data Compared to 3pt Moving Average Filter')
xlim([0,(length(p)-1)/300])
xlabel('Time (sec)')
ylabel('Pulsatile Blood Flow')
legend('Raw','Filtered')
grid
hold off

figure
plot(t,p, '--')
hold on
plot(t,pfiltered)
xlim([.1,.35])
ylim([8, 10])
title({'Effect of Filter on Raw Data'; '(First Peak Zoomed In)'})
xlabel('Time (sec)')
ylabel('Pulsatile Blood Flow')
legend('Raw','Filtered')
grid
hold off

clear all

%Question 5f

H = ones(31);
for m = 1:length(H(:,1))
    H(m,:) = 2^(-(m-1));
end

for i = 2:length(H(:,1)) %recycled code from #3
    Hi = H(:,i-1);
    Hi = circshift(Hi,1);
    Hi(1) = 0;
    H(:,i) = Hi;
end

t = linspace(0,30,31);
impres = zeros(1,31);
impres(1,1) = 1;
pimp = H*impres';

un = ones(1,31);
pun = H*un';

figure
plot(t,pimp)
hold on 
plot(t,pun)

ax = gca;
xlim([-5,30])
xlabel('Time (days)')
ax.XAxisLocation = 'origin'; 
ylim([-1,3])
ylabel('Protein Concentration')
ax.YAxisLocation = 'origin';
title({'Protein Concentration'; '(With and Without Protein Synthesis)'})
legend('Degradation only','Degradation and Synthesis','Location','south')
grid
hold off