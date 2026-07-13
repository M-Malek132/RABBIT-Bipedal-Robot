function [ N ] = BSpline( n,p,u )
% n = number of datapoints
% p = degree
% u = corresponding control point

%initialize

N = zeros(1,n+1); % coefs of Cp
m = n+p+1; % Count of knot vectors

% generating knot vectors
U =zeros(1,m+1);
U(:,p+1:n+2) = linspace(0,1,m+1-2*p );
U(:,n+2:m+1) = 1;

%speciall cases
if u== U(1)
    N(1) = 1.0;
    return;

else  if u==U(m+1)
        N(n+1) = 1.0;

        return;
end
end
%degree of 0 coef
for k=1:m+1
    if U(k) > u
        N(k-1) = 1.0;
        k = k-1;
        break;
    else if U(k) == u
            N(k) = 1.0;
            break;
    end
    end
end
% k=k;

%main loop
for d=1:p % finding coefs in desired degree

    N(k-d) = ((U(k+1) - u) / (U(k+1) - U(k-d+1))) * N(k-d+1);% calculate south west corner term only
    for i= k-d+1:k-1
        N(i) = ((u-U(i)) / (U(i+d) - U(i))) * N(i) + ((U(i+d+1) -u) /(U(i+d+1) - U(i+1)))*N(i+1) ;
    end
    N(k) = ((u-U(k))/( U(k+d) - U(k)  ))*N(k); % calculate north west corner term only
end
end