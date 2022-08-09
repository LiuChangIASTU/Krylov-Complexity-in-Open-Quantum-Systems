% % % ---------------------------------------------------------------------
% % % READ ME.
 
% % % This code evaluate the spectrum of the nonHermitian hopping model and
% % % plots the FIG.4 in the main text.
 
% % % (1), ns are the same definition in the main text.
 
% % % (2), nmax is the numerical truncation of the half-infinite chain. Here we
% % % choose nmax=2000 and ns=1000 to make sure that nmax is large enough so
% % % that it donot affect the low energy physics.
 
% % % (3), alpha_b, beta_b, alpha_d beta_d, are the same definition as in the
% % % main text, and their value are the same as those in caption of FIG.4.
 
% % % (4), gamma is the same definition as in main text. For blue dashed line,
% % % we take gamma=0, for red solid line, we take gamma=0.1 which is in the
% % % gapped phase.

% % % ---------------------------------------------------------------------
% % % This part is some settings for plotting a figure, such as the position of 
% % % subplot, the fontsize and so on, which is irrelevant of physics.
figure('Color','White');len=0.328;labelsize=28;axissize=16;
xdis=[0.18 0.1 0.18 0.1]-0.05;
ydis=[0.12 0.12 0.17 0.17]-0.05;
subpo=[0.83,0.64]-0.05;

annotation('textbox',[0.045 0.85 0.1 0.1], 'String','${\rm (a)}$','EdgeColor','none','FontSize',24,'Interpreter','Latex');
annotation('textbox',[0.51-0.03 0.85 0.1 0.1], 'String','${\rm (b)}$','EdgeColor','none','FontSize',24,'Interpreter','Latex');
annotation('textbox',[0.045 0.4 0.1 0.1], 'String','${\rm (c)}$','EdgeColor','none','FontSize',24,'Interpreter','Latex');
annotation('textbox',[0.51-0.03 0.4 0.1 0.1], 'String','${\rm (d)}$','EdgeColor','none','FontSize',24,'Interpreter','Latex');
% % % ---------------------------------------------------------------------
% % % This part plot FIG.4(a), and the definition for the parameter in below
% % % are the same as those in caption of FIG.4, Here we choose gamma=0.01
% % % which is in the gapless phase.
nmax=2000;ns=1000;
alpha_b=0.3435;beta_b=0.66287;
alpha_d=0.35189;beta_d=2.8126;
gamma=0.007;
% % % the below are detailed code for plot.

alpha_a=gamma*alpha_d;beta_a=gamma*beta_d;
L=L_gen(alpha_a,beta_a,alpha_b,beta_b,nmax,ns);
[V1,D1]=eig(1i*L);
eigval1=diag(D1);
[~,ind_sort]=sort(real(eigval1),'descend');
Vplt1=zeros(nmax+1,nmax+1);
eigvalplt1=zeros(nmax+1,1);
for n=1:(nmax+1)
    Vplt1(:,n)=V1(:,ind_sort(n));
    eigvalplt1(n)=eigval1(ind_sort(n));
end
clear V1;clear W1;clear eigval1;clear D1;clear temp

axes('position',[0+xdis(1),0.5+ydis(1),len,len])
plot(imag(eigvalplt1),-real(eigvalplt1)+real(eigvalplt1(1)),'bo');
set(gca,'TickLabelInterpreter','Latex','FontSize',axissize);
    xlabel('$\varepsilon''$','Interpreter','Latex','FontSize',labelsize);
    ylabel('$\varepsilon''''$','Interpreter','Latex','FontSize',labelsize);
xlim([-20,20]);ylim([0,0.3]);

% % % ---------------------------------------------------------------------
% % % This part plot FIG.4(b), and the definition for the parameter in below
% % % are the same as those in caption of FIG.4, Here we choose gamma=0.1
% % % which is in the gapped phase.
nmax=2000;ns=1000;
alpha_b=0.3435;beta_b=0.66287;
alpha_d=0.35189;beta_d=2.8126;
gamma=0.04;
% % % the below are detailed code for plot.

alpha_a=gamma*alpha_d;beta_a=gamma*beta_d;
L=L_gen(alpha_a,beta_a,alpha_b,beta_b,nmax,ns);

[V1,D1]=eig(1i*L);
eigval1=diag(D1);
[~,ind_sort]=sort(real(eigval1),'descend');

Vplt2=zeros(nmax+1,nmax+1);
eigvalplt2=zeros(nmax+1,1);

for n=1:(nmax+1)
    Vplt2(:,n)=V1(:,ind_sort(n));
    eigvalplt2(n)=eigval1(ind_sort(n));
end
clear V1;clear W1;clear eigval1;clear D1;clear temp

axes('position',[0.5+xdis(2),0.5+ydis(2),len,len])
plot(imag(eigvalplt2),-real(eigvalplt2)+real(eigvalplt2(1)),'bo');
set(gca,'TickLabelInterpreter','Latex','FontSize',axissize);
xlabel('$\varepsilon''$','Interpreter','Latex','FontSize',labelsize);
ylabel('$\varepsilon''''$','Interpreter','Latex','FontSize',labelsize);
xlim([-15,15]);ylim([0,2.2]);
    
% % % ---------------------------------------------------------------------
% % % This part plot FIG.4(c), and the definition for the parameter in below
% % % are the same as FIG.4(b)

axes('position',[0+xdis(3),0+ydis(3),len,len])
index=1;plot(0:nmax,reshape(abs(Vplt2(:,index)).^2,[1,nmax+1]),'b-o','MarkerSize',5,'DisplayName',['$','|\varphi_{',num2str(index),'}|^{2}','$']);hold on;
index=2;plot(0:nmax,reshape(abs(Vplt2(:,index)).^2,[1,nmax+1]),'r-','LineWidth',1.5,'DisplayName',['$','|\varphi_{',num2str(index),'}|^{2}','$']);
set(gca,'TickLabelInterpreter','Latex','FontSize',axissize);
xlabel('$n$','Interpreter','Latex','FontSize',labelsize);
ylabel('$|\phi(n)|^{2}$','Interpreter','Latex','FontSize',22);
ylim([0,0.009]);    
% % % ---------------------------------------------------------------------
% % % This part plot FIG.4(d) and the process is more complicated compared
% % % with those above. Given ns, to find the corresponding gamma_c, we need to
% % % do procedure described in below:
 
% % % Fix ns, change nmax, find gamma_c for each nmax, and read out the
% % % saturation value of gamma_c when nmax goes to infinity. 
 
% % % In practice, we find that gamma_c already saturates when nmax=2*ns and remains 
% % % the same for larger nmax. 
 
% % % The value of ns we are stored in variable nlst in below, and the
% % % corresponding saturated gamma_c's value are stored in variable data. The
% % % linear regression result of log(gamma_c) versus log(ns) are stored in
% % % variable temp4, and reader can read out the power-law exponent: -0.8648
% % % used in the main text from temp4(1);

nslst=[300,500,750,1000,1250,1500,2000];
gammacdata=[0.06364,0.04167,0.03,0.02333,0.01889,0.01611,0.01278];

temp4=regress(transpose(log(gammacdata)),[transpose(log(nslst)),ones(length(nslst),1)]);
axes('position',[0.5+xdis(4),0+ydis(4),len,len])
loglog(nslst,exp(temp4(1)*log(nslst)+temp4(2)*ones(1,length(nslst))),'r-','LineWidth',1.5);hold on;
loglog(nslst,gammacdata,'bo','MarkerSize',14);hold on;
set(gca,'TickLabelInterpreter','Latex','FontSize',axissize);
xlabel('$n_{s}$','Interpreter','Latex','FontSize',labelsize);
ylabel('$\gamma_{c}$','Interpreter','Latex','FontSize',labelsize);

% % % ---------------------------------------------------------------------
% % % This part is an example :Fix ns, get different gamma_c for different nmax

% % % alpha_b=0.3435;beta_b=0.66287;
% % % alpha_d=0.35189;beta_d=2.8126;
% % 
% % % nslst=[500,750,1000,1250,1500,2000];Nns=length(nslst);
% % % gammalst=linspace(0,0.055,100);Ngamma=length(gammalst);
% % % Nnmax=2;
% % % DDgap=zeros(Nns,Nnmax,Ngamma-2);
% % % gapp=zeros(Nns,Nnmax,Ngamma);
% % % 
% % % for ins=1:Nns
% % % ns=nslst(ins);
% % % nmaxlst=[2*ns+200,2*ns+500];
% % % for inmax=1:Nnmax
% % % 
% % % nmax=nmaxlst(inmax);
% % %     for igamma=1:Ngamma
% % %         tic;['ns=',num2str(ns),', inmax=',num2str(inmax),', igamma=',num2str(igamma)]
% % %         gamma=gammalst(igamma);
% % %         alpha_a=gamma*alpha_d;beta_a=gamma*beta_d;
% % %         L=L_gen(alpha_a,beta_a,alpha_b,beta_b,nmax,ns);
% % %         eigval1=eig(1i*L);
% % %         [eigval,ind_sort]=sort(real(eigval1),'descend');
% % %         gapp(ins,inmax,igamma)=real(eigval(2)-eigval(3));
% % %         toc
% % %     end
% % %         
% % %         temp=(gapp(ins,inmax,1:(Ngamma-2))+gapp(ins,inmax,3:Ngamma)-2*gapp(ins,inmax,2:(Ngamma-1)))/gammalst(2);
% % %         DDgap(ins,inmax,:)=abs(temp);
% % % 
% % % end
% % % end
% % 
% % % ins=6;
% % % ns=nslst(ins);
% % % nmaxlst=[2*ns+200,2*ns+500];
% % % for inmax=1:Nnmax
% % %     figure('Color','White');
% % %     subplot(1,2,1)
% % %     plot(gammalst(2:(Ngamma-1)),reshape(DDgap(ins,inmax,:),[1,Ngamma-2]),'b-o','DisplayName',['$n_{max}=',num2str(nmaxlst(inmax)),', n_{s}=',num2str(ns),'$']);hold on; 
% % %     legend('Interpreter','Latex');
% % %     subplot(1,2,2)
% % %     plot(gammalst,reshape(gapp(ins,inmax,:),[1,Ngamma]),'r-o')
% % % end
% % % ---------------------------------------------------------------------

% % %----------------------------------------------------------------------

% % % -------This function generate the hopping matrix with nearest
% neighbor hopping's bn and onsite damping.
function L=L_gen(alpha_a,beta_a,alpha_b,beta_b,nmax,ns)
L=zeros(nmax+1,nmax+1);
L(1,1)=1i*beta_a;
for n=1:ns
    an=1i*(alpha_a*n+beta_a);
    bn=alpha_b*n+beta_b;
    L(n+1,n+1)=an;
    L(n,n+1)=bn;L(n+1,n)=bn;
end
if ns<nmax
for n=(ns+1):nmax

    L(n+1,n+1)=an;
    L(n,n+1)=bn;L(n+1,n)=bn;
end
end

end
