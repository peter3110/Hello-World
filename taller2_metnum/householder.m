function [X] = householder(A,b)
%HOUSEHOLDER Summary of this function goes here
%   Detailed explanation goes here
   n = size(A, 1);
   R = A;
   Q_tr = eye(size(A));
   for i = 1:n-1
      x = [zeros(i-1,1); R(i:n,i)];

      e_i = zeros(n, 1);
      e_i(i) = 1;

      alpha = norm(x);
      
      e_i = e_i*alpha;
      
      u = x - e_i;

      v = u/norm(u);

      R = R - 2*v*(v'*R);
      Q_tr = Q_tr - 2*v*(v'*Q_tr) ;
   end

   Q = Q_tr';
   X = R\(Q'*b);
   return
end
