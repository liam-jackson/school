load acetylene
plotmatrix([x1 x2])
X = [x1 x2];
D = x2fx(X,'interaction');
D(:,1) = []; % No constant term
k = 0:1e-5:5e-3;
B = ridge(y,D,k);
figure
plot(k,B,'LineWidth',2)
ylim([-100 100])
grid on 
xlabel('Ridge Parameter') 
ylabel('Standardized Coefficient') 
title('Ridge Trace') 
legend('x1','x2','x3','x1x2','x1x3','x2x3')
