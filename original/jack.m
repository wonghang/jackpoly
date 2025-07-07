% function f=jack(l,x,alpha). Computes the Schur function corresponding
% to partition l with variables specified in x.
%
% Plamen Koev, 2003

function f=jack(l,xx,alpha1);
  global x ja Lp Lmax n lma alpha
  x=xx;
  alpha=alpha1;
  lma=[];
  
  n=length(x);

  Lp=0;
  while (length(l)>Lp)&(l(Lp+1)>0)
     Lp=Lp+1;
  end % number of parts in l
  l=l(1:Lp);
  Lmax=l+1;

  ja(prod(Lmax)+1,n)=inf; % allocate space
  if isequal(class(xx(1)),'sym')
    ja=sym(ja);
  end
 
  if Lp>0
     lma(Lp)=1;
     for i=Lp-1:-1:1
        lma(i)=lma(i+1)*Lmax(i+1);
     end

     initialize(1,0);
     f=jack1(n,0,nmu(l),nmu(l));
  else
     f=1;
  end
  
  clear x Lmax ja Lp n lma jlambda alpha;

function f=initialize(k,l);
  global Lmax Lp ja n lma

  if k<=Lp
     m=Lmax(k)-1;
     if (k>1)
        m=min(m,part(l,k-1));
     end
     for i=1:m
        l=l+lma(k);
        ja(l+1,1:n)=inf*ones(1,n);
        initialize(k+1,l);
     end
  end

% given an integer pn that represents a partition, this function returns
% part i of the partition in O(1) time

function f=part(pn,i)
global lma

if i>length(lma)
   f=0;
else
   if i==1
      f=floor(pn/lma(i));
   else
      f=floor(mod(pn,lma(i-1))/lma(i));
   end
end


% nmu computes the unique integer that represents the partition l

function f=nmu(l)
global Lmax Lp
  f=0;
  for i=1:Lp
     f=Lmax(i)*f;
     if i<=length(l)
        f=f+l(i);
     end
  end


% Given vector x and partition, represented by an integer l, jack1(j,k,lambda,l)
% computes the Jack polynomial J_lambda^alpha (x_1,...,x_j), where the partitions
% l, such that lambda-l is horizontal strip, are generated to have 
% lambda(i)=l(i) for i=1,2,...,k. The parameter k is designed to keep 
% track of the recursion. from outside the function should be called
% with ja(j,0,lambda,lambda). This function can only be called after an initialization
% by the function jack

function f=jack1(j,k,lambda,l);
  global x Lp Lmax n ja lma alpha

  s=1;

  if (1<=j)&(l>0)
    t=ja(l+1,j);
    if (k==0)&(t~=inf)
      s=t;
    elseif part(l,j+1)>0
      s=0;
    elseif j==1
      s=x(1)^part(l,1)*prod(1+alpha*[0:part(l,1)-1]);
    else
      i=k+(k==0);
      s=jack1(j-1,0,l,l)*AB(lambda,l)*x(j)^lm(lambda,l); 
      while part(l,i)>0
        if part(l,i)>part(l,i+1)
          if part(l,i)>1
            s=s+jack1(j,i,lambda,l-lma(i));
          else
            s=s+jack1(j-1,0,l-lma(i),l-lma(i))*AB(lambda,l-lma(i))*x(j)^lm(lambda,l-lma(i));
          end
        end
        i=i+1;
      end
    end
    if (k==0)
      ja(l+1,j)=s;
    end
  end 
  f=s;

function f=AB(l,m)
global Lp alpha
  f=1;
  for i=1:Lp
     for j=1:part(m,i)
        if l_t(l,j)==l_t(m,j)
           f=f/(l_t(m,j)-i+alpha*(part(m,i)-j+1));
        else
           f=f/(l_t(m,j)-i+1+alpha*(part(m,i)-j));
        end
     end
  end
  
  for i=1:Lp
     for j=1:part(l,i)
        if l_t(l,j)==l_t(m,j)
           f=f*(l_t(l,j)-i+alpha*(part(l,i)-j+1));
        else
           f=f*(l_t(l,j)-i+1+alpha*(part(l,i)-j));
        end
     end
  end


% computes |lambda/mu|
function f=lm(l,m)
global Lp
  f=0;
  for i=1:Lp
     f=f+part(l,i)-part(m,i);
  end
  
   
% compute the q-th part of the transpose of a partition lt(l,q)=#(l(i)|l(i)>=q)
function f=l_t(l,q)

i=1;
f=0;
while (part(l,i)>=q) 
   f=f+1;
   i=i+1;
end

