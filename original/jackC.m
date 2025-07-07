% function [j,Rx]=jackC(lambda,x,alpha)
% computes the C normalization of the Jack polynomial 
%
% Rx is a vector of all Jack polynomials for 
% partitions mu \subset lambda that are computed in the
% process anyway.
%
% by Plamen Koev, December 16, 2016

function [j,Rx]=jackC(lam,x,alpha)
  
  MAX=sum(lam);
  p=[];
  q=[];
  if length(lam)==length(x)
      x=[x 0];
  end
  

  % ss are the partial sums
  
  n=length(x);
  lambda=lam;
    
  Lp=length(lambda);
  while lambda(Lp)==0
     Lp=Lp-1;
  end % number of parts in l
  lambda=lambda(1:Lp);
  
  % figure out number of partitions lambda with |lambda|<=MAX
  % and number of parts less than or equal to Lp
  f=1:MAX+1; 
  for i=2:Lp-1
      for j=i+1:MAX+1
          f(j)=f(j)+f(j-i);
      end
  end
  f=f(end); % only need to allocate space for partitions of <=n-1 parts
  XY=0;
%
%  Let N(kappa) be the index for partition kappa
%  D(N(kappa)) stores N((kappa,1))
%  Sx(i,1:n) stores various jack functions of X
%  with the first entry S_kappa(x_1), second entry
%  S_kappa(x_1,x_2), and so on.
%  Sy is similar to Sx but it is for Y.
% 
  D=zeros(f,1);
  Sx=zeros(f,n);
  Sx(1,:)=1;
  Rx=Sx(:,n);
% xn(i,j)=x(i)^(j-1), yn(i,j)=y(i)^(j-1)
  xn=ones(n,MAX+1);
  if size(x,2)>1, x=x.'; end
  for i=2:MAX+1
     xn(:,i)=xn(:,i-1).*x;
  end
  prodx=cumprod(x);
  
  if XY
     Sy=Sx;
     yn=ones(n,MAX+1);
     if size(y,2)>1, y=y.'; end; 
     for i=2:MAX+1
         yn(:,i)=yn(:,i-1).*y;
     end
     prody=cumprod(y);
  end
  
  l=zeros(1,Lp);
  z=ones(1,Lp);
  kt=-(1:Lp);   % kt(i)=alpha*l(i)-i, it is updated instead of recomputed every time
  
  cc1=0;
  ss=zeros(1,MAX+1); % these are the partial sums
  ss(1)=1;
  
  sl=1;  % sl=sum(l), it is updated instead of recomputed every time
%
% Note that h is the number of nonzero parts of l,
% w is N(l), i.e., the index for partition l.
% heap tells us where next block of memory is available after the current
% l(h)=1,..,m is stored, where m is the number of possible elements
% for l(h).
%
  h=1;
  ww=ones(1,Lp);
  heap=lambda(1)+2;  
  d = zeros(1,Lp);   % d=l-[l(2:Lp) 0], it is the places we will be adding boxes
  while h>0
     if (l(h)<lambda(h))&&(h==1||l(h)<l(h-1))&&(MAX>=sl)&&(z(h)~=0)
        l(h)=l(h)+1;
        if l(h)==1&&h>1&&h<n
           D(ww(h))=heap;
           ww(h)=heap;
           m=min(lambda(h),MAX-sl+l(h));
           heap=heap+min(m,l(h-1));   
        else
           ww(h)=ww(h)+1;
        end
        w=ww(h);
%   Update Q
        c=-(h-1)/alpha+l(h)-1;
        zn=prod(p+c)*alpha;
        dn=prod(q+c)*(kt(h)+h+1);
        if XY
           zn=zn*alpha*l(h);
           dn=dn*(n+alpha*c);
           for j=1:h-1
               delta=kt(j)-kt(h);
               zn=zn*delta;
               dn=dn*(delta-1);
           end
        end
        kt(h)=kt(h)+alpha;
        for j=1:h-1
            delta=kt(j)-kt(h);
            zn=zn*delta;
            dn=dn*(delta+1);
        end
        z(h)=z(h)*zn/dn;
        sl=sl+1;        
%
%    Only need to do work when number of nonzero
%    parts in l is less than n.
%
        if h<n
           if h>1
              d(h-1)=d(h-1)-1;
           end
           d(h)=l(h);
           cc=prod(h+1-alpha+kt(1:h))/prod(h+kt(1:h));
           % pp is the index of l-ones(1,Lp)
           pp=l(1);
           k=2;
           while k<=h&&l(k)>1
              pp=D(pp)+l(k)-2;
              k=k+1;
           end
           Sx(w,h)=cc*prodx(h)*Sx(pp,h);
           if XY        
              Sy(w,h)=cc*prody(h)*Sy(pp,h);
           end
           g=find(d(1:h)>0);   % g are the indices of places we will be adding 
           lg=length(g);
           slm=1;              % this is sum(l-mu)+1 which we will be updating
           nhstrip=prod(d(g)+1)-1;  % number of mu that are horizontal strips of l
           mu=l;
           mt=kt;              % mt(i)=alpha*m(i)-i   
           blm=ones(1,lg);     % This is a vector to store the beta_{lambda,mu} update
           lmd=l(g)-d(g);      % This gives us the lower bound of mu
           for i=1:nhstrip
               j=lg;
               gz=g(lg);    
               while mu(gz)==lmd(j)
                  mu(gz)=l(gz);
                  mt(gz)=kt(gz);
                  slm=slm-d(gz);
                  j=j-1;
                  gz=g(j);  
               end
  
               t=kt(gz)-mt(gz);
               blm(j)=blm(j)*(1+t);
               dn=(t+alpha);
               for r=1:gz-1
                   q1=mt(r)-mt(gz);
                   q2=kt(r)-mt(gz);
                   blm(j)=blm(j)*(q1+alpha-1)*(1+q2);
                   dn=dn*q1*(alpha+q2);
               end
               blm(j)=blm(j)/dn;
               mu(gz)=mu(gz)-1;
               mt(gz)=mt(gz)-alpha;
               slm=slm+1;
%            
%    Push blm(j) all the way to the end because subsequent partitions
%    of mu will be smaller than the current one. 
%  
               if j<lg      
                  blm(j+1:end)=blm(j);
               end
%
%    Find out the index for partition mu. 
% 
               nmu=mu(1)+1;
               for k=2:h-(mu(h)==0)
                   nmu=D(nmu)+mu(k)-1;
               end
            
               for k=h+1:n
                   Sx(w,k)=Sx(w,k)+blm(j)*Sx(nmu,k-1)*xn(k,slm);
               end
               if XY
                  for k=h+1:n
                      Sy(w,k)=Sy(w,k)+blm(j)*Sy(nmu,k-1)*yn(k,slm);
                  end
               end
           end
           for k=h+1:n
               Sx(w,k)=Sx(w,k)+Sx(w,k-1);
           end
          
              ss(sl)=ss(sl)+z(h)*Sx(w,n);
              Rx(w)=z(h)*Sx(w,n) * ...
                      hooklengths(l,alpha,3) ...
                      /alpha^(sum(l));
        end
        if h<Lp,         % increase the length of the partition if possible
           z(h+1)=z(h);  % carry over Q_kappa and index pointer
           ww(h+1)=w;
           h=h+1;
        end
     else
        sl=sl-l(h);
        l(h)=0;
        kt(h)=-h;
        h=h-1;
     end
  end
  
  % recover jack_lam
  
  pp=lam(1)+1;
  k=2;
  while k<=length(lam)
     pp=D(pp)+lam(k)-1;
     k=k+1;
  end
  
  j=Rx(pp); 
    % multiplying by the upper hook lengths to get the J normalization
