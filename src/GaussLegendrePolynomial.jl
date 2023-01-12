module GaussLegendrePolynomial

export GauLeg!
export GauLeg
export GauLeg_Log
export GauLeg_Log!
export Legfit_array_cmplx!
export Legev_cmplx!

EPS=1e-11

"""
Gauss-Legendre polynomical. The return value is stored in x and w.
"""
function GauLeg!(::Type{T},a::T,b::T,x::Vector{T},w::Vector{T},eps=EPS)where{T<:AbstractFloat}
  n=length(x)
  m=(n+1)÷2
  xm=0.5*(a+b); xl=0.5*(b-a)
  for i =0:m-1
    z=cos(π*(i+0.75)/(n+0.5))
    z1 = 1e100
    pp=0.0
    while abs(z-z1)>eps
      p1=1.0; p2=0.0
      for j=0:n-1
        p3=p2; p2=p1
        p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1)
      end
      pp=n*(z*p1-p2)/(z*z-1.0)
      z1=z
      z=z1-p1/pp
    end
    x[i+1]=xm-xl*z  
    x[n-i]=xm+xl*z
    w[i+1]=2.0*xl/((1.0-z*z)*pp*pp)   
    w[n-i]=w[i+1]
  end 
  return nothing
end # function GauLeg!

function GauLeg!(lo,hi,x,w,eps=EPS)
  GauLeg!(typeof(lo),lo,hi,x,w,eps)
end

function GauLeg_Log!(lo,hi,x,w,eps=EPS)
  GauLeg!(log10(lo), log10(hi), x,w);
  for i=eachindex(x)
    x[i]=10^(x[i])
    w[i]=x[i]*w[i]*log(10.0)
  end
end

function GauLeg(::Type{T},a,b,n::Integer,eps=EPS)where{T<:AbstractFloat}
  x=zeros(T,n); w=zeros(T,n)

  m=(n+1)÷2
  xm=0.5*(a+b); xl=0.5*(b-a)
  for i =0:m-1
    z=cos(π*(i+0.75)/(n+0.5))
    z1 = 1e100
    pp=0.0
    while abs(z-z1)>eps
      p1=1.0; p2=0.0
      for j=0:n-1
        p3=p2; p2=p1
        p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1)
      end
      pp=n*(z*p1-p2)/(z*z-1.0)
      z1=z
      z=z1-p1/pp
    end
    x[i+1]=xm-xl*z  
    x[n-i]=xm+xl*z
    w[i+1]=2.0*xl/((1.0-z*z)*pp*pp)   
    w[n-i]=w[i+1]
  end 
  return x,w
end # function GauLeg

function GauLeg(lo,hi,n::Integer,eps=EPS)
  x,w=GauLeg(typeof(lo),lo,hi,n)
  return x,w
end

function GauLeg_Log(lo,hi,n::Integer,eps=EPS)
  x,w=GauLeg(log10(lo), log10(hi), n);
  for i=1:n
    x[i]=10^(x[i])
    w[i]=x[i]*w[i]*log(10.0)
  end
  return x,w
end


"""
Poly_P() returns the value of l-order Legendre polynomical at point x
(l阶勒让德多项式在x点的值)
"""
function Poly_P(l::Integer,x)     
  P0=1;
  P1=x;
  P2=0;

  for i=1:l-1
      P2=((2*i+1)*x*P1-i*P0)/(i+1);
      P0=P1;
      P1=P2;
	end


  if l==0
      return P0;
  elseif l==1
      return P1;
  else
      return P2;
  end
end

"""
对函数func()在区间(a,b)上的一系列点func[]，这些点取值在相应多项式零点上，对其做勒让德逼近，返回多项式展开的系数
"""
function Legfit_array_cmplx!(func::Vector{Complex{T}},c::Vector{Complex{T}},n::Integer)where{T<:AbstractFloat}
  x_lo::T = -1.0
  x_hi::T = 1.0
	x,w = GauLeg(x_lo, x_hi, n)
	
	for j=1:n
	  sum::Complex{T}=0.0;
	  for k=1:n
	    sum+=w[k]*func[k]*Poly_P(j-1,x[k]);
    end
	  c[j]=sum*(2*(j-1)+1)/2.0;
	end


end

function Legev_cmplx!(a::T, b::T, c::Vector{Complex{T}},m::Integer, x::T)where{T<:AbstractFloat}
	sum::Complex{T} =0.0
	xr=(b-a)/2.0
	xm=(b+a)/2.0
	
	for i=1:m
	  sum += c[i]*Poly_P(i-1,(x-xm)/xr)
  end
	return sum
end


end
