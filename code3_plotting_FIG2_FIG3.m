% % % ---------------------------------------------------------------------
% % % READ ME
% % % This code plots FIG.2 and FIG.3 in the main text, before plotting, one
% % % shall import the data from fermion model and SYK model, since the data
% % % are calculated in different matlab-code. 
load('syk_12_data.mat');
load('fermion_13_data.mat');
% % % ---------------------------------------------------------------------
% % % This part plot FIG.2

figure('Color','White')
len1=0.32;len2=0.77;fontsize=26;

axes('position',[0.17,0.15,len1,len2])
imagesc(D(:,:,2),'CdataMapping','scaled');
xlabel('$n$','Interpreter','Latex','FontSize',fontsize);
ylabel('$m$','Interpreter','Latex','FontSize',fontsize);
set(gca,'TickLabelInterpreter','Latex','FontSize',fontsize);
colorbar('northoutside');

axes('position',[0.55,0.15,len1,len2])
imagesc(Dmn,'CdataMapping','scaled');
xlabel('$n$','Interpreter','Latex','FontSize',fontsize);
set(gca,'TickLabelInterpreter','Latex','FontSize',fontsize);
colorbar('northoutside');

annotation('textbox',[0.1 0.85 0.1 0.1], 'String','${\rm (a)}$','EdgeColor','none','FontSize',27,'Interpreter','Latex');
annotation('textbox',[0.49 0.85 0.1 0.1], 'String','${\rm (b)}$','EdgeColor','none','FontSize',27,'Interpreter','Latex');

% % % ---------------------------------------------------------------------
% % % This part plot FIG.3(a)

figure('Color','White')
len1=0.67;len2=0.37;fontsize=20;
annotation('textbox',[0.02 0.90 0.1 0.1], 'String','${\rm (a)}$','EdgeColor','none','FontSize',23,'Interpreter','Latex');
annotation('textbox',[0.02 0.41 0.1 0.1], 'String','${\rm (b)}$','EdgeColor','none','FontSize',23,'Interpreter','Latex');

axes('position',[0.15,0.60,len1,len2])
plot(0:nmax,real(2*N-Dnn(:,1)),'g-o','DisplayName',['$O=','\psi_{1}',', M=','\psi_{i}','$']);hold on;
plot(0:nmax,real(Dnn(:,2)),'r-*','DisplayName',['$O=','i\psi_{1}\psi_{2}',', M=','\psi_{i}','$']);hold on;
plot(0:nmax,real(Dnn(:,3)),'b-^','DisplayName',['$O=','\psi_{1}',', M=','i\psi_{i}\psi_{j}','$']);hold on;
plot(0:nmax,real(Dnn(:,4)),'c-h','DisplayName',['$O=','i\psi_{1}\psi_{2}',', M=','i\psi_{i}\psi_{j}','$']);hold on;
text(7,19,'$\uparrow$','Interpreter','Latex','FontSize',28);
text(7,16,'${\rm n_{s}}$','Interpreter','Latex','FontSize',25);
legend('Interpreter','Latex','FontSize',15);

set(gca,'TickLabelInterpreter','Latex','FontSize',20);
ylabel('${\rm d_{n}}$','Interpreter','Latex','FontSize',25);


% % % ---------------------------------------------------------------------
% % % This part plot FIG3.(b)
axes('position',[0.15,0.13,len1,len2])

plot(1:50,ff,'g-o','MarkerSize',5,'DisplayName',['$','O=','\psi_{1}',', M=','\psi_{i}','$']);hold on;
plot(1:50,bf,'r-*','MarkerSize',5,'DisplayName',['$','O=','n_{1}',', M=','\psi_{i}','$']);hold on;
plot(1:50,fb,'b-^','MarkerSize',5,'DisplayName',['$','O=','\psi_{1}',', M=','2n_{i}','$']);hold on;
plot(1:50,bb,'c-h','MarkerSize',5,'DisplayName',['$','O=','n_{1}',', M=','2n_{i}','$']);hold on;
set(gca,'TickLabelInterpreter','Latex','FontSize',20);
text(24,6,'$\uparrow$','Interpreter','Latex','FontSize',28);
text(24,4.7,'${\rm n_{s}}$','Interpreter','Latex','FontSize',24);
legend('Interpreter','Latex','FontSize',14);
xlabel('${\rm n}$','Interpreter','Latex','FontSize',27);
ylabel('${\rm d_{n}}$','Interpreter','Latex','FontSize',25);
% % % ---------------------------------------------------------------------
