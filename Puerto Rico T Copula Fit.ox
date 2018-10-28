#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.oxh>

// PDF of the GPD distribution
densgpd(x,xi,sigma,u){
return((1/sigma)*(1+(xi/sigma)*(x-u))^(-(1/xi+1)));
}

// CDF of the GPD distribution
probgpd(x,xi,sigma,u){
 return(1-(1+(xi/sigma)*(x-u))^(-1/xi));
}

// Posterior distribution of the T-copula
postt(data, rho, p,v, xi1, sigma1, u1, mu1, eta1, p1, xi2, sigma2, u2, mu2, eta2, p2,  a1, b1, c1, d1, meanu1, sdu1, a2, b2, c2, d2, meanu2, sdu2){
	 decl cont, score, g1,g2;
     score=0;
     for (cont=0;cont<rows(data);cont++){	 
		if (data[cont][0]<= u1 && data[cont][1] <= u2){
			g1 = quant(sumr(p1.*probgamma(data[cont][0],eta1,eta1./mu1)),v);  g2 = quant(sumr(p2.*probgamma(data[cont][1],eta2,eta2./mu2)),v);
			score += log((1/sqrt(1-rho^2))*(gammafact((v+2)/2)*gammafact(v/2)/gammafact((v+1)/2)^2)*(1+g1^2/v+g2^2/v+(g1^2*g2^2)/(v^2))^((v+1)/2)/(1+(g1^2+g2^2-2*rho*g1*g2)/(v*(1-rho^2)))^((v+2)/2)) + log(sumr(p1.*densgamma(data[cont][0],eta1,eta1./mu1)))+ log(sumr(p2.*densgamma(data[cont][1],eta2,eta2./mu2)));	}
		else if (data[cont][0]<= u1 && data[cont][1] > u2){
			g1 = quant(sumr(p1.*probgamma(data[cont][0],eta1,eta1./mu1)),v); g2 = quant(sumr(p2.*probgamma(u2,eta2,eta2./mu2))+(1-sumr(p2.*probgamma(u2,eta2,eta2./mu2)))*probgpd(data[cont][1],xi2,sigma2,u2),v);
	     score += log((1/sqrt(1-rho^2))*(gammafact((v+2)/2)*gammafact(v/2)/gammafact((v+1)/2)^2)*(1+g1^2/v+g2^2/v+(g1^2*g2^2)/(v^2))^((v+1)/2)/(1+(g1^2+g2^2-2*rho*g1*g2)/(v*(1-rho^2)))^((v+2)/2))  + log(sumr(p1.*densgamma(data[cont][0],eta1,eta1./mu1))) +log(1-sumr(p2.*probgamma(u2,eta2,eta2./mu2)))+log(densgpd(data[cont][1],xi2,sigma2,u2));
			}
		else if (data[cont][0]> u1 && data[cont][1] <= u2){
			g2 = quant(sumr(p2.*probgamma(data[cont][1],eta2,eta2./mu2)),v); g1 = quant(sumr(p1.*probgamma(u1,eta1,eta1./mu1))+(1-sumr(p1.*probgamma(u1,eta1,eta1./mu1)))*probgpd(data[cont][0],xi1,sigma1,u1),v);
			score += log((1/sqrt(1-rho^2))*(gammafact((v+2)/2)*gammafact(v/2)/gammafact((v+1)/2)^2)*(1+g1^2/v+g2^2/v+(g1^2*g2^2)/(v^2))^((v+1)/2)/(1+(g1^2+g2^2-2*rho*g1*g2)/(v*(1-rho^2)))^((v+2)/2))  +log(1-sumr(p1.*probgamma(u1,eta1,eta1./mu1)))+log(densgpd(data[cont][0],xi1,sigma1,u1)) + log(sumr(p2.*densgamma(data[cont][1],eta2,eta2./mu2)));   
			}					
		else if (data[cont][0]> u1 && data[cont][1] > u2){
			g2 = quant(sumr(p2.*probgamma(u2,eta2,eta2./mu2))+(1-sumr(p2.*probgamma(u2,eta2,eta2./mu2)))*probgpd(data[cont][1],xi2,sigma2,u2),v); g1 = quant(sumr(p1.*probgamma(u1,eta1,eta1./mu1))+(1-sumr(p1.*probgamma(u1,eta1,eta1./mu1)))*probgpd(data[cont][0],xi1,sigma1,u1),v);
			score += log((1/sqrt(1-rho^2))*(gammafact((v+2)/2)*gammafact(v/2)/gammafact((v+1)/2)^2)*(1+g1^2/v+g2^2/v+(g1^2*g2^2)/(v^2))^((v+1)/2)/(1+(g1^2+g2^2-2*rho*g1*g2)/(v*(1-rho^2)))^((v+2)/2))  +log(1-sumr(p1.*probgamma(u1,eta1,eta1./mu1)))+log(densgpd(data[cont][0],xi1,sigma1,u1)) + log(1-sumr(p2.*probgamma(u2,eta2,eta2./mu2)))+log(densgpd(data[cont][1],xi2,sigma2,u2));	   
			}
		}
	score += sumr((c1-1).*log(eta1)-d1.*eta1-(b1+1).*log(mu1)-a1./mu1)+log((sigma1^(-1)*(1+xi1)^(-1)*(1+2*xi1)^(-0.5)))-0.5*((u1-meanu1)/sdu1)^2+sumr((c2-1).*log(eta2)-d2.*eta2-(b2+1).*log(mu2)-a2./mu2)+log((sigma2^(-1)*(1+xi2)^(-1)*(1+2*xi2)^(-0.5)))-0.5*((u2-meanu2)/sdu2)^2;	
    score += sumr(0.5.*log(v./(v+3))+0.5.*log(polygamma(v./2,1)-polygamma((v+1)./2,1)-2.*(v+3)/(v.*(v+1).^2)));
	return(score);	
	}

	main(){
//MCMC algorithm choices	
decl it = 25000, burn=5000, thin=20, num=(it-burn)/thin;
ranseed(11);

// Starting values of the chain
decl  mu2 = <100,300>, eta2 = <2,2>, xi2 = 0.5, u2 = 463, sigma2 = 10, p2 = <0.5,0.5>,  dim2 = columns(mu2);
decl mu1 = <10,200>, eta1 = <2,2>, xi1 = 0.1, u1 = 373, sigma1 = 10, p1 = <0.5,0.5>,  dim1 = columns(mu1);
decl p = <1>, rho = <0.9>, v= <2>, dim = columns(rho);

// Load the dataset - Use the right directory
decl data = loadmat("/Users/manuele/Dropbox/Dani/Data_River/Puerto Rico Fit Dataset.dat");

// Values for the prior distribution
decl b1 = constant(2.1,1,dim1), a1 = constant(meanr(mu1),1,dim1).*(b1[0]-1), d1 = constant(0.5,1,dim1),  c1 = d1[0]./constant(meanr(eta1),1,dim1), meanu1 =373, sdu1 = 50;
decl b2 = constant(2.1,1,dim2), a2 = constant(meanr(mu2),1,dim2).*(b2[0]-1), d2 = constant(0.5,1,dim2),  c2 = d2[0]./constant(meanr(eta2),1,dim2), meanu2 = 463, sdu2 = 50;

// Creation of the chains
decl ximc1 = zeros(it+1,1), sigmamc1 = zeros(it+1,1), umc1 = zeros(it+1,1), mumc1 = zeros(it+1,dim1), etamc1 = zeros(it+1,dim1), pmc1=zeros(it+1,dim1);
decl ximc2 = zeros(it+1,1), sigmamc2 = zeros(it+1,1), umc2 = zeros(it+1,1), mumc2 = zeros(it+1,dim2), etamc2 = zeros(it+1,dim2), pmc2 = zeros(it+1,dim2);
decl pmc = zeros(it+1,dim), rhomc = zeros(it+1,dim), vmc= zeros(it+1,1);

// Variables counting the numbers of accepted steps
decl accxi1 = 0, accsigma1 = 0, accu1 = 0, accmu1 = constant(0,1,dim1), acceta1 = constant(0,1,dim1);
decl accxi2 = 0, accsigma2 = 0, accu2 = 0,  accmu2 = constant(0,1,dim2), acceta2 = constant(0,1,dim2);
decl accrho = constant(0,1,dim), accv=constant(0,1,1);

// Variances of the proposal distributions
decl Vxi1 = 0.25, Vsigma1 = 5, Vu1 = 10, Vp1 = 50, Vmu1 = <5,5>, Veta1 = <1,1>;
decl Vxi2 = 0.25, Vsigma2 = 5, Vu2 = 10, Vp2 = 50, Vmu2 = <5,5>, Veta2 = <1,1>;
decl Vp = 50, Vrho = constant(0.05,1,dim), Vv=<0.5>;

// Definition of variables for the algorithm
decl M=maxc(data), vp, vprob, vaux, restrict, xip1, xip2, sigmap1, sigmap2, up1, up2, pp, pp1, pp2,  xiprob1, xiprob2, sigmaprob1, unifo, sigmaprob2, uprob1, uprob2, pprob, pprob1, pprob2, adj, i, m, Ga, Ga1, Ga2, j, mup1, muaux1, mup2, muaux2, etap1, etaaux1, etap2, etaaux2, gammaprob1, gammaprob2, restrict1, restrict2, rhop, rhoaux, rhoprob, muprob1, muprob2, etaprob1,etaprob2;

// Initialization of the chains
ximc1[0]=xi1; sigmamc1[0]=sigma1; umc1[0]=u1; mumc1[0][]=mu1; etamc1[0][]=eta1; pmc1[0][]=p1;
ximc2[0]=xi2; sigmamc2[0]=sigma2; umc2[0]=u2; mumc2[0][]=mu2; etamc2[0][]=eta2; pmc2[0][]=p2;
pmc[0][]=p; rhomc[0][]=rho;	 vmc[0]=v;

						  
 for (i= 1; i <= it; ++i){

// Simulation of rho
	 rhop=rhomc[i-1]+Vrho*rann(1,1);
	 while(rhop>1 || rhop< -1){rhop=rhomc[i-1]+Vrho*rann(1,1);}
	 rhoprob=exp(postt(data,rhop,pmc[i-1][],vmc[i-1],ximc1[i-1],sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2)-postt(data,rhomc[i-1],pmc[i-1][],vmc[i-1], ximc1[i-1],sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2));
	 if (ranu(1,1)<rhoprob){rhomc[i]=rhop; accrho+=1;} else {rhomc[i]=rhomc[i-1]; }

// Simulation of p (only one copula component)
	 pmc[i][]=pmc[i-1][];

// Simulation of v
	 vp= max(rangamma(1,1,vmc[i-1]^2/Vv^2,vmc[i-1]/Vv^2),0.001);
	 vprob=exp(postt(data,rhomc[i][],pmc[i][],vp,ximc1[i-1],sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2)-postt(data,rhomc[i][],pmc[i][],vmc[i-1], ximc1[i-1],sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2)- log(densgamma(vp,vmc[i-1].^2 ./Vv.^2,vmc[i-1]./Vv.^2))+log(densgamma(vmc[i-1],vp.^2 ./Vv.^2,vp./Vv.^2)));
	 if(ranu(1,1)<vprob) {vmc[i]=vp; accv+=1;} else {vmc[i]=vmc[i-1];}
	
// Simulation of xi1
 	 xip1=max(ximc1[i-1]+Vxi1*rann(1,1),-4.4999);
	 while (xip1<(-(sigmamc1[i-1])/(M[0]-umc1[i-1]))) {xip1=max(ximc1[i-1]+Vxi1*rann(1,1),-4.4999);}
	 xiprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i],xip1,sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i-1],sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )+log(probn((ximc1[i-1]+sigmamc1[i-1]/(M[0]-umc1[i-1]))/Vxi1))-log(probn((xip1+sigmamc1[i-1]/(M[0]-umc1[i-1]))/Vxi1)));
	 if(ranu(1,1)<xiprob1) {ximc1[i]=xip1; accxi1+=1;} else {ximc1[i]=ximc1[i-1];}
 		
// Simulation of sigma1
 	 if(ximc1[i]<0){
		sigmap1=max(sigmamc1[i-1]+Vsigma1*rann(1,1),0.0001);
		while(sigmap1< (-ximc1[i]*(M[0]-umc1[i-1]))){sigmap1=max(sigmamc1[i-1]+Vsigma1*rann(1,1),0.0001); }
		sigmaprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmap1,umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) + log(probn((sigmamc1[i-1]+ximc1[i]*(M[0]-umc1[i-1]))/Vsigma1))-log(probn((sigmap1+ximc1[i]*(M[0]-umc1[i-1]))/Vsigma1)));}
	 else{
		sigmap1=max(rangamma(1,1,sigmamc1[i-1]^2/Vsigma1^2,sigmamc1[i-1]/Vsigma1^2),0.0001);
		sigmaprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmap1,umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i-1],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) + log(densgamma(sigmamc1[i-1],sigmap1^2/Vsigma1^2,sigmap1/Vsigma1^2))-log(densgamma(sigmap1,sigmamc1[i-1]^2/Vsigma1^2,sigmamc1[i-1]/Vsigma1^2)));}
	 if(ranu(1,1)<sigmaprob1){sigmamc1[i]=sigmap1; accsigma1+=1;} else {sigmamc1[i]=sigmamc1[i-1];}		

// Simulation of u1
	 if(ximc1[i]<0){m=M[0]+sigmamc1[i]/ximc1[i];} else {m=minc(data)[0];}
	 up1=umc1[i-1]+Vu1*rann(1,1);
	 while (up1<m) {up1=umc1[i-1]+Vu1*rann(1,1);}
	 uprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],up1,mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i-1],mumc1[i-1][],etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) + log(probn((umc1[i-1]-m)/Vu1))-log(probn((up1-m)/Vu1)));
	 if(ranu(1,1)<uprob1){umc1[i]=up1; accu1+=1;} else {umc1[i]=umc1[i-1];}

// Simulation of mu1 and eta1	
 	 mup1=mumc1[i-1][];muaux1=mumc1[i-1][];	 	
	 mup1[0]=max(rangamma(1,1,mumc1[i-1][0]^2/Vmu1[0]^2,mumc1[i-1][0]/Vmu1[0]^2),0.001);
	 muprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mup1,etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],muaux1,etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(mup1[0],muaux1[0]^2/Vmu1[0]^2,muaux1[0]/Vmu1[0]^2))+log(densgamma(muaux1[0],mup1[0]^2/Vmu1[0]^2,mup1[0]/Vmu1[0]^2)));
	 if(ranu(1,1)<muprob1) {muaux1[0]=mup1[0];} else {mup1[0]=muaux1[0];}
     mup1[1]=max(rangamma(1,1,mumc1[i-1][1]^2/Vmu1[1]^2,mumc1[i-1][1]/Vmu1[1]^2),0.001);
     muprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mup1,etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],muaux1,etamc1[i-1][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(mup1[1],muaux1[1]^2/Vmu1[1]^2,muaux1[1]/Vmu1[1]^2))+log(densgamma(muaux1[1],mup1[1]^2/Vmu1[1]^2,mup1[1]/Vmu1[1]^2)));
     if (ranu(1,1)<muprob1) {muaux1[1]=mup1[1];} else {mup1[1]=muaux1[1];}
     mumc1[i][]=mup1;
     restrict1=(mumc1[i][]-sortr(mumc1[i][]))*(mumc1[i][]-sortr(mumc1[i][]))';
     if(restrict1>0) {mumc1[i][]=mumc1[i-1][];}
     for (j=0; j<dim1;++j) { if (mumc1[i][j]-mumc1[i-1][j]!=0){accmu1[j]+=1;}}

	 etap1=etamc1[i-1][]; etaaux1=etamc1[i-1][];
	 etap1[0]=max(rangamma(1,1,etamc1[i-1][0]^2/Veta1[0]^2,etamc1[i-1][0]/Veta1[0]^2),0.001);
	 etaprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etap1,pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etaaux1,pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(etap1[0],etaaux1[0]^2/Veta1[0]^2,etaaux1[0]/Veta1[0]^2))+log(densgamma(etaaux1[0],etap1[0]^2/Veta1[0]^2,etap1[0]/Veta1[0]^2)));
	 if (ranu(1,1)<etaprob1) {etaaux1[0]=etap1[0]; acceta1[0]+=1;} else { etap1[0]=etaaux1[0];}
	 etap1[1]=max(rangamma(1,1,etamc1[i-1][1]^2/Veta1[1]^2,etamc1[i-1][1]/Veta1[1]^2),0.001);
	 etaprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etap1,pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etaaux1,pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(etap1[1],etaaux1[1]^2/Veta1[1]^2,etaaux1[1]/Veta1[1]^2))+log(densgamma(etaaux1[1],etap1[1]^2/Veta1[1]^2,etap1[1]/Veta1[1]^2)));
     if(ranu(1,1)<etaprob1) {etaaux1[1]=etap1[1]; acceta1[1]+=1;} else { etap1[1]=etaaux1[1];}
     etamc1[i][]=etap1;

// Simulation of p1
	 Ga1=zeros(1,dim1);
	 for(j=0;j<dim1;j++) {Ga1[j]=max(rangamma(1,1,Vp1.*pmc1[i-1][j],1),0.001);}
		pp1=Ga1./sumr(Ga1);
		pprob1=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pp1,ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i-1][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) +sumr(log(pmc1[i-1][]))-sumr(log(pp1)));
	 if(ranu(1,1)<pprob1){pmc1[i][]=pp1;} else {pmc1[i][]=pmc1[i-1][];}
	
// Simulation of xi2
	 xip2=max(ximc2[i-1]+Vxi2*rann(1,1),-4.4999);
	 while (xip2<(-(sigmamc2[i-1])/(M[1]-umc2[i-1]))) {xip2=max(ximc2[i-1]+Vxi2*rann(1,1),-4.4999);}
	 xiprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],xip2,sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i-1],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )+log(probn((ximc2[i-1]+sigmamc2[i-1]/(M[1]-umc2[i-1]))/Vxi2))-log(probn((xip2+sigmamc2[i-1]/(M[1]-umc2[i-1]))/Vxi2)));
	 if(ranu(1,1)<xiprob2) {ximc2[i]=xip2; accxi2+=1;} else {ximc2[i]=ximc2[i-1];}
		
// Simulation of sigma2
	 if(ximc2[i]<0){
		sigmap2=max(sigmamc2[i-1]+Vsigma2*rann(1,1),0.0001);
		while(sigmap2< (-ximc2[i]*(M[1]-umc2[i-1]))){sigmap2=max(sigmamc2[i-1]+Vsigma2*rann(1,1),0.0001); }
		sigmaprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmap2,umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) + log(probn((sigmamc2[i-1]+ximc2[i]*(M[1]-umc2[i-1]))/Vsigma2))-log(probn((sigmap2+ximc2[i]*(M[1]-umc2[i-1]))/Vsigma2)));}
	 else{
		sigmap2=max(rangamma(1,1,sigmamc2[i-1]^2/Vsigma2^2,sigmamc2[i-1]/Vsigma2^2),0.0001);
		sigmaprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmap2,umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i-1],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) + log(densgamma(sigmamc2[i-1],sigmap2^2/Vsigma2^2,sigmap2/Vsigma2^2))-log(densgamma(sigmap2,sigmamc2[i-1]^2/Vsigma2^2,sigmamc2[i-1]/Vsigma2^2)));}
	 if(ranu(1,1)<sigmaprob2){sigmamc2[i]=sigmap2; accsigma2+=1;} else {sigmamc2[i]=sigmamc2[i-1];}		
	  	
// Simulation of u
	 if(ximc2[i]<0){m=M[1]+sigmamc2[i]/ximc2[i];} else {m=minc(data)[1];}
	 up2=umc2[i-1]+Vu2*rann(1,1);
	 while (up2<m) {up2=umc2[i-1]+Vu2*rann(1,1); }
	 uprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],up2,mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i-1],mumc2[i-1][],etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) + log(probn((umc2[i-1]-m)/Vu2))-log(probn((up2-m)/Vu2)));
	 if(ranu(1,1)<uprob2){umc2[i]=up2; accu2+=1;} else {umc2[i]=umc2[i-1];}

// Simulation of mu2 and eta2
	 mup2=mumc2[i-1][];muaux2=mumc2[i-1][];	 	
	 mup2[0]=max(rangamma(1,1,mumc2[i-1][0]^2/Vmu2[0]^2,mumc2[i-1][0]/Vmu2[0]^2),0.001);
	 muprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mup2,etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],muaux2,etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(mup2[0],muaux2[0]^2/Vmu2[0]^2,muaux2[0]/Vmu2[0]^2))+log(densgamma(muaux2[0],mup2[0]^2/Vmu2[0]^2,mup2[0]/Vmu2[0]^2)));
	 if(ranu(1,1)<muprob2) {muaux2[0]=mup2[0];} else {mup2[0]=muaux2[0];}
	 mup2[1]=max(rangamma(1,1,mumc2[i-1][1]^2/Vmu2[1]^2,mumc2[i-1][1]/Vmu2[1]^2),0.001);
	 muprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mup2,etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],muaux2,etamc2[i-1][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(mup2[1],muaux2[1]^2/Vmu2[1]^2,muaux2[1]/Vmu2[1]^2))+log(densgamma(muaux2[1],mup2[1]^2/Vmu2[1]^2,mup2[1]/Vmu2[1]^2)));
     if(ranu(1,1)<muprob2) {muaux2[1]=mup2[1];} else {mup2[1]=muaux2[1];}
	 mumc2[i][]=mup2;
	 restrict2=(mumc2[i][]-sortr(mumc2[i][]))*(mumc2[i][]-sortr(mumc2[i][]))';
	 if(restrict2>0) {mumc2[i][]=mumc2[i-1][];}
     for(j=0; j<dim2;++j) { if (mumc2[i][j]-mumc2[i-1][j]!=0){accmu2[j]+=1;}}

	 etap2=etamc2[i-1][]; etaaux2=etamc2[i-1][];
	 etap2[0]=max(rangamma(1,1,etamc2[i-1][0]^2/Veta2[0]^2,etamc2[i-1][0]/Veta2[0]^2),0.001);
	 etaprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mumc2[i][],etap2,pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mumc2[i][],etaaux2,pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(etap2[0],etaaux2[0]^2/Veta2[0]^2,etaaux2[0]/Veta2[0]^2))+log(densgamma(etaaux2[0],etap2[0]^2/Veta2[0]^2,etap2[0]/Veta2[0]^2)));
	 if (ranu(1,1)<etaprob2) {etaaux2[0]=etap2[0]; acceta2[0]+=1;} else { etap2[0]=etaaux2[0];}
	 etap2[1]=max(rangamma(1,1,etamc2[i-1][1]^2/Veta2[1]^2,etamc2[i-1][1]/Veta2[1]^2),0.001);
	 etaprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mumc2[i][],etap2,pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mumc2[i][],etaaux2,pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )- log(densgamma(etap2[1],etaaux2[1]^2/Veta2[1]^2,etaaux2[1]/Veta2[1]^2))+log(densgamma(etaaux2[1],etap2[1]^2/Veta2[1]^2,etap2[1]/Veta2[1]^2)));
	 if(ranu(1,1)<etaprob2) {etaaux2[1]=etap2[1]; acceta2[1]+=1;} else { etap2[1]=etaaux2[1];}
	 etamc2[i][]=etap2;

// Simulation of p2
     Ga2=zeros(1,dim2);
	 for(j=0;j<dim2;j++) {Ga2[j]=max(rangamma(1,1,Vp2.*pmc2[i-1][j],1),0.001);}
		pp2=Ga2./sumr(Ga2);
		pprob2=exp(postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mumc2[i][],etamc2[i][],pp2,a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 )-postt(data,rhomc[i][],pmc[i][],vmc[i], ximc1[i],sigmamc1[i],umc1[i],mumc1[i][],etamc1[i][],pmc1[i][],ximc2[i],sigmamc2[i],umc2[i],mumc2[i][],etamc2[i][],pmc2[i-1][],a1,b1,c1,d1,meanu1,sdu1,a2,b2,c2,d2,meanu2,sdu2 ) +sumr(log(pmc2[i-1][]))-sumr(log(pp2)));
	 if(ranu(1,1)<pprob2){pmc2[i][]=pp2;} else {pmc2[i][]=pmc2[i-1][];}
	
	
// Adjustment of Proposal Variances
     if(imod(i,50)==0){
        adj=min(0.01,(trunc(i/50)+1)^(-0.5));
        if (accxi1>22){Vxi1=Vxi1*exp(adj);} else {Vxi1=Vxi1*exp(-adj);}
		if (accsigma1>22){Vsigma1=Vsigma1*exp(adj);} else {Vsigma1=Vsigma1*exp(-adj);}
		if (accu1>22){Vu1=Vu1*exp(adj);} else {Vu1=Vu1*exp(-adj);}
		if (accxi2>22){Vxi2=Vxi2*exp(adj);} else {Vxi2=Vxi2*exp(-adj);}
		if (accsigma2>22){Vsigma2=Vsigma2*exp(adj);} else {Vsigma2=Vsigma2*exp(-adj);}
		if (accu2>22){Vu2=Vu2*exp(adj);} else {Vu2=Vu2*exp(-adj);}
		if (accv>22){Vv=Vv*exp(adj);} else {Vv=Vv*exp(-adj);}
		for (j=0; j<dim; ++j){ 	  if (accrho[j]>22){Vrho[j]=Vrho[j]*exp(adj);} else {Vrho[j]=Vrho[j]*exp(-adj);} accrho[j]=0;}
		for (j=0; j<dim2; ++j){
		    if (accmu2[j]>22){Vmu2[j]=Vmu2[j]*exp(adj);} else {Vmu2[j]=Vmu2[j]*exp(-adj);}
			 if (acceta2[j]>22){Veta2[j]=Veta2[j]*exp(adj);} else {Veta2[j]=Veta2[j]*exp(-adj);}
			accmu2[j]=0; acceta2[j]=0;}
		    if (accmu1[0]>22){Vmu1[0]=Vmu1[0]*exp(adj);} else {Vmu1[0]=Vmu1[0]*exp(-adj);}
			if (acceta1[0]>22){Veta1[0]=Veta1[0]*exp(adj);} else {Veta1[0]=Veta1[0]*exp(-adj);}
			accmu1[0]=0; acceta1[0]=0; accxi1=0; accsigma1=0; accu1=0; accxi2=0; accsigma2=0; accu2=0; accv=0;}			

// Print iteration
	 println(i,  " ", vmc[i], " ", accv, " ", Vv);
 }

// Thinning of the chain  
decl gpd1=zeros(num,3), gammamu1=zeros(num,2), gammaeta1=zeros(num,2), gammap1=zeros(num,2), gpd2=zeros(num,3), gammamu2=zeros(num,2), gammaeta2=zeros(num,2), gammap2=zeros(num,2), dep=zeros(num,2),  k;
for (k=0; k<num; ++k){
gpd1[k][0]=ximc1[burn+(k+1)*thin], gpd1[k][1]=sigmamc1[burn+(k+1)*thin], gpd1[k][2]=umc1[burn+(k+1)*thin];
gammamu1[k][0]=mumc1[burn+(k+1)*thin][0], gammamu1[k][1]=mumc1[burn+(k+1)*thin][1] ;
gammaeta1[k][0]=etamc1[burn+(k+1)*thin][0], gammaeta1[k][1]=etamc1[burn+(k+1)*thin][1] ;
gammap1[k][0]=pmc1[burn+(k+1)*thin][0], gammap1[k][1]=pmc1[burn+(k+1)*thin][1] ;
gpd2[k][0]=ximc2[burn+(k+1)*thin], gpd2[k][1]=sigmamc2[burn+(k+1)*thin], gpd2[k][2]=umc2[burn+(k+1)*thin];
gammamu2[k][0]=mumc2[burn+(k+1)*thin][0], gammamu2[k][1]=mumc2[burn+(k+1)*thin][1] ;
gammaeta2[k][0]=etamc2[burn+(k+1)*thin][0], gammaeta2[k][1]=etamc2[burn+(k+1)*thin][1] ;
gammap2[k][0]=pmc2[burn+(k+1)*thin][0], gammap2[k][1]=pmc2[burn+(k+1)*thin][1] ;
dep[k][0]=rhomc[burn+(k+1)*thin][0],dep[k][1]=vmc[burn+(k+1)*thin][0];
}
savemat("GPD1.mat",gpd1);
savemat("GPD2.mat",gpd2);
savemat("MU1.mat",gammamu1);
savemat("ETA1.mat",gammaeta1);
savemat("P1.mat",gammap1);
savemat("MU2.mat",gammamu1);
savemat("ETA2.mat",gammaeta1);
savemat("P2.mat",gammap1);
savemat("DEP.mat",dep);						  

}
