function [ne,curp] = from_particles(ne,curp,ni,x1,dx,px,pv,nm,pweight)

% Accumulate particle density and proton current density to grid.

ne(:) = 1.e-20;
curp(:) = 1.e-20;

for m= 1:nm

% ri is where the particle is on the grid.  ri is a real 
%  number, and ii and ip are the integer grid positions to the left and 
%  right of the particle location.

  ri= ( px(m,2)-x1 )/dx + 1;
  ii= floor(ri);
  ip= ii + 1;

% The particle pweight contributes to density at both points ii and ip
%  with weights fii and fip.  This is called linear weighting.  Note 
%  that if ri=ii (particle is at lower ii grid point), fii=1 and 
%  fip=0, while if the particle is at ri=ip, fip=1 and fii=0.

  fip= ri - ii;
  fii= 1. - fip;
  ne(ii)= ne(ii) + fii*pweight;
  ne(ip)= ne(ip) + fip*pweight;
  curp(ii,:)= curp(ii,:) + fii*pweight*pv(m,:,2);
  curp(ip,:)= curp(ip,:) + fip*pweight*pv(m,:,2);

end;  % end of 'for m='

% Fix density overflow into buffer region.

ne(   2)= ne(   2) + ne(ni); 
ne(ni-1)= ne(ni-1) + ne( 1);
curp(   2,:)= curp(   2,:) + curp(ni,:);
curp(ni-1,:)= curp(ni-1,:) + curp( 1,:);

% Fix buffer region assuming periodic boundary conditions.

ne( 1)= ne(ni-1);
ne(ni)= ne(   2);
curp( 1,:)= curp(ni-1,:);
curp(ni,:)= curp(   2,:);

end
