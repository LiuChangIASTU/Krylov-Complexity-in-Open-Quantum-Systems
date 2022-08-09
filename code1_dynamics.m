% % % ---------------------------------------------------------------------
% % % READ ME.

% % % This code evaluate the dynamics of K(t), and plots the FIG.1 in the main
% % % text. We first explain the variables in the code:
 
% % % (1), ns are the same definition in the main text.
 
% % % (2), nmax is the numerical truncation of the half-infinite chain. Here we
% % % choose nmax=1000 and ns=300 to make sure that nmax is large enough so
% % % that such numerical truncation donot affect the dynamics of K(t).
 
% % % (3), alpha_b, beta_b, alpha_d, beta_d are the same definition as in the
% % % main text.
 
% % % (4), gamma is the same definition as in main text. For blue dashed line,
% % % we take gamma=0, for red solid line, we take gamma=0.1 which is in the
% % % gapped phase.

% % % ---------------------------------------------------------------------
% % % Definition of the variables mentioned in the READ ME. Readers are welcome
% % % to try other values.

nmax=1000;ns=300;

alpha_b=0.4162;beta_b=3.295;
alpha_d=0.2067;beta_d=2.851;
gammalst=[0,0.1];

% % % ---------------------------------------------------------------------

phi0=[1;zeros(nmax,1)];

T1=6;T2=T1;
tlst1=linspace(0,T1,2000);tlst2=linspace(T1,T2,0);
Nt=length(tlst1)+length(tlst2);
tlst=[tlst1,tlst2];

Ngamma=length(gammalst);

Krylov=zeros(Ngamma,Nt);

igamma=1;gamma=gammalst(igamma);
alpha_a=gamma*alpha_d;beta_a=gamma*beta_d;
L=L_gen(alpha_a,beta_a,alpha_b,beta_b,nmax,ns);
Krylov(igamma,:)=Kflow_gen(phi0,L,tlst1,tlst2);

igamma=2;gamma=gammalst(igamma);
alpha_a=gamma*alpha_d;beta_a=gamma*beta_d;
L=L_gen(alpha_a,beta_a,alpha_b,beta_b,nmax,ns);
Krylov(igamma,:)=Kflow_gen(phi0,L,tlst1,tlst2);


figure('Color','White');fontsize=20;

igamma=1;gamma=gammalst(igamma);
plot(tlst,reshape(Krylov(igamma,:),[1,Nt]),'b--','LineWidth',1.7,'DisplayName',['$','\gamma=',num2str(gamma),'$']);hold on
igamma=2;gamma=gammalst(igamma);
plot(tlst,reshape(Krylov(igamma,:),[1,Nt]),'r-','LineWidth',1.7,'DisplayName',['$','\gamma=',num2str(gamma),'$']);hold on

set(gca,'TickLabelInterpreter','Latex');
xlabel('$t$','Interpreter','Latex','FontSize',fontsize);
ylabel('$\mathcal{K}(t)$','Interpreter','Latex','FontSize',fontsize);
ylim([0,500]);
xlim([0,5.4]);
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

% % % -------This function evaluate the Krylov complexity given the wave
% function on the half-chain.
function y=Kcomp(phi)
phi=abs(phi).^2;
phi=phi/sum(phi);
y=sum(transpose(0:(length(phi)-1)).*phi);


end

% % % -------This function refresh the small timestep of the wavefunction
% by using 4-th-order's RungeKutta method
function y=RungeKutta(yn,L,deltat)
k1=deltat * 1i*L*(yn);
k2=deltat * 1i*L*(yn+0.5*k1);
k3=deltat * 1i*L*(yn+0.5*k2);
k4=deltat * 1i*L*(yn+k3);
y=yn+k1/6+k2/3+k3/3+k4/6;
end

% % % -------This function generate the full Krylov complexity dependence
% on a time series.
function y=Kflow_gen(phi0,L,tlst1,tlst2)
Nt1=length(tlst1);
deltat1=tlst1(2)-tlst1(1);
K1=zeros(Nt1,1);
phi=phi0;
for it1=1:Nt1
%     tic;[num2str(it1)]
        phi=RungeKutta(phi,L,deltat1);
        K1(it1)=Kcomp(phi);
%     toc
end

Nt2=length(tlst2);
K2=zeros(Nt2,1);
[V,D,W]=eig(L);
temp=(W'*phi0)./diag(W'*V);eigval=diag(D);
for it=1:Nt2
%     tic;it
        t=tlst2(it);
        phit=sum(V*diag(exp(1i*t*eigval).*temp),2);
        K2(it)=Kcomp(phit);
%     toc
end
y=[K1;K2];
end

