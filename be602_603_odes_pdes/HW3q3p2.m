%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Defining Meshgrid for part 2
x2 = -1:.1:1;
y2 = x2; 
z2 = x2;

[X2, Y2, Z2] = meshgrid(x2,y2,z2);

R2 = sqrt(X2.^2 + Y2.^2 + Z2.^2);
theta2 = acos(Z2./R2);
phi2 = atan2(Y2,X2);
xdumb2 = cos(theta2); 
xdumb2sq = 1-xdumb2.^2; 

%Dirichlet Radius of Sphere
b2 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining arrays for part 3 
x3 = -1:.01:1;

phi3 = 0;
xdumb3 = x3;
xdumb3sq = 1-xdumb3.^2;

%Dirichlet Radius of Sphere
b3 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining Meshgrid for part 5
x5 = -5:0.1:5;
y5 = x5; 
z5 = x5;

[X5, Y5, Z5] = meshgrid(x5,y5,z5);

R5 = sqrt(X5.^2 + Y5.^2 + Z5.^2);
theta5 = acos(Z5./R5);
phi5 = atan2(Y5,X5);
xdumb5 = cos(theta5);
xdumb5sq = 1-xdumb5.^2;

%Dirichlet Radius of Sphere
b5 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Final Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:1:length(Lallow)
%     for j = 1:1:length(Mallow)
%         L = Lallow(i);
%         M = Mallow(j);
%         PvalForSum = Pval(i,j);
%         PvalForSum
%         PvalForSum2(i,j,:,:,:) = subs(PvalForSum,xdumb,xdumb2);
%         PvalForSum3(i,j,:) = subs(PvalForSum,xdumb,xdumb3);
%         PvalForSum5(i,j,:,:,:) = subs(PvalForSum,xdumb,xdumb5);
%     end
% end

%Here we go:

%Summation for Final Solution
% for CLindex = 1:1:length(Lallow)
%     for CMindex = 1:1:length(Mallow)
%         L = Lallow(CLindex);
%         M = Mallow(CMindex);
%         if M <= L
%             T2 = T2 ...
%                 + (c(CLindex,CMindex)...
%                 .*((b2./R2).^(L+1))... 
%                 .*squeeze(PvalForSum2(CLindex,CMindex,:,:,:))...
%                 .*cos(M.*phi2));
%             T3 = T3 ...
%                 + (c(CLindex,CMindex)...
%                 .*((1./b3).^(L+1))... 
%                 .*squeeze(PvalForSum3(CLindex,CMindex,:))...
%                 .*cos(M.*phi3));
%             T5 = T5 ...
%                 + (c(CLindex,CMindex)...
%                 .*((b5./R5).^(L+1))... 
%                 .*squeeze(PvalForSum5(CLindex,CMindex,:,:,:))...
%                 .*cos(M.*phi5));
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simplify(Pval);
%part 2 soln
tic 
s2_00 = c(1,1).*((b2./R2).^(0+1)).*1;
s2_20 = c(3,1).*((b2./R2).^(2+1)).*((1/2).*(3.*xdumb2.^2-1));
s2_40 = c(5,1).*((b2./R2).^(4+1)).*((1/8).*(35.*xdumb2.^4-30.*xdumb2.^2+3)); 
s2_60 = c(7,1).*((b2./R2).^(6+1)).*((1/16).*(231.*xdumb2.^6-315.*xdumb2.^4+105.*xdumb2.^2-5)); 
s2_22 = c(3,3).*((b2./R2).^(2+1)).*((3.*xdumb2sq));
s2_42 = c(5,3).*((b2./R2).^(4+1)).*(15/2.*xdumb2sq.*((7.*xdumb2.^2)-1)); 
s2_62 = c(7,3).*((b2./R2).^(6+1)).*(105/8.*xdumb2sq.*(33*xdumb2.^4-18*xdumb2.^2+1)); 

m0terms2 = s2_00 + s2_20 + s2_40 + s2_60; 
m2terms2 = s2_22 + s2_42 + s2_62; 

T2 = m0terms2 + (m2terms2 .* cos(2*phi2));
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part 3 soln
tic
s3_00 = c(1,1).*((1./b3).^(0+1)).* 1;
s3_20 = c(3,1).*((1./b3).^(2+1)).* (1/2).*(3.*xdumb3.^2-1);
s3_40 = c(5,1).*((1./b3).^(4+1)).* (1/8).*(35.* xdumb3.^4-30.*xdumb3.^2+3); 
s3_60 = c(7,1).*((1./b3).^(6+1)).* (1/16).*(231.*xdumb3.^6-315.*xdumb3.^4+105.*xdumb3.^2-5); 
s3_22 = c(3,3).*((1./b3).^(2+1)).* (3.*xdumb3sq);
s3_42 = c(5,3).*((1./b3).^(4+1)).* (15/2.*xdumb3sq.*((7.*xdumb3.^2)-1)); 
s3_62 = c(7,3).*((1./b3).^(6+1)).* 105/8.*xdumb3sq.*(33*xdumb3.^4-18*xdumb3.^2+1); 

m0terms3 = s3_00 + s3_20 + s3_40 + s3_60; 
m2terms3 = s3_22 + s3_42 + s3_62; 

T3 = m0terms3 + (m2terms3 .* cos(2*phi3));
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part 5 soln
tic
s5_00 = c(1,1).*((b5./R5).^(0+1)).* 1;
s5_20 = c(3,1).*((b5./R5).^(2+1)).* (1/2).*(3.*xdumb5.^2-1);
s5_40 = c(5,1).*((b5./R5).^(4+1)).* (1/8).*(35.* xdumb5.^4-30.*xdumb5.^2+3); 
s5_60 = c(7,1).*((b5./R5).^(6+1)).* (1/16).*(231.*xdumb5.^6-315.*xdumb5.^4+105.*xdumb5.^2-5); 
s5_22 = c(3,3).*((b5./R5).^(2+1)).* (3.*xdumb5sq);
s5_42 = c(5,3).*((b5./R5).^(4+1)).* (15/2.*xdumb5sq.*((7.*xdumb5.^2)-1)); 
s5_62 = c(7,3).*((b5./R5).^(6+1)).* 105/8.*xdumb5sq.*(33*xdumb5.^4-18*xdumb5.^2+1); 

m0terms5 = s5_00 + s5_20 + s5_40 + s5_60; 
m2terms5 = s5_22 + s5_42 + s5_62; 

T5 = m0terms5 + (m2terms5 .* cos(2*phi5));
toc
