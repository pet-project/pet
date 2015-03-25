function [cur,ve,e] = from_fields(cur,ve,e,b,bx0,curp,ne,ni,nim,dx2)

% Calculate current
for i= 2:ni-1
  cur(i,2) = -( b(i+1,3,2) - b(i-1,3,2) )/dx2;
  cur(i,3) =  ( b(i+1,2,2) - b(i-1,2,2) )/dx2;
end
cur(ni,2)= cur(  2,2); cur(ni,3)= cur(  2,3);
cur( 1,2)= cur(nim,2); cur( 1,3)= cur(nim,3);

% Calculate ve
for iv= 1:3
  ve(:,iv) = ( curp(:,iv) - cur(:,iv) ) ./ ne(:);
end

% Calculate E
e(:,1) = - ve(:,2).*b(:,3,2) + ve(:,3).*b(:,2,2);
e(:,2) = - ve(:,3).*bx0      + ve(:,1).*b(:,3,2);
e(:,3) = - ve(:,1).*b(:,2,2) + ve(:,2).*bx0     ;

end
