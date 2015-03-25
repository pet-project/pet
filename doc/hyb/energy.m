function [ben,pen,ten] = energy(it,b,ni,nid,by0,pv,nm,pweight,pen1,i_subtract_thermal_energy )

% Calculate b field energy ben (per grid point).

  ben= 0.;
  for i= 2:ni-1
    ben= ben + (b(i,2,2)-by0)^2 + b(i,3,2)^2;
  end;
  ben= .5*ben/nid;      % Normalized to get average (perturbed) B energy density.
        
% Calculate particle kinetic energy

  pen= 0;
  for m= 1:nm
    pen= pen + pv(m,1,2)^2 + pv(m,2,2)^2 + pv(m,3,2)^2;
  end;
  pen= 0.5*pweight*pen/nid;
  if( i_subtract_thermal_energy==1 & it>1 )
    pen= pen - pen1;
  end
  ten= ben + pen;
  
end