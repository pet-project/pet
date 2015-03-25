function [px,pv] = step_particles_boris( px,pv,nm,e,b,ni,nim,lsim,x1,dx,bx0,dt,dtd2,iboris,ipc )

pe= zeros(3,1);
pb= zeros(3,1);
if( iboris==1 )
  pvm= zeros(3,1);
  pvs= zeros(3,1);
  pvp= zeros(3,1);
  t= zeros(3,1);
  s= zeros(3,1);
end

for m= 1:nm
  
  ri= ( px(m)-x1 )/dx + 1;
  ii= floor(ri);
  ip= ii + 1;
  fip= ri - ii;
  fii= 1. - fip;

% Now get e and b field felt by particle, pe.

  pe= fii*e(ii,:) + fip*e(ip,:);
  pb(2)= fii*b(ii,2) + fip*b(ip,2);
  pb(3)= fii*b(ii,3) + fip*b(ip,3);
  if( iboris==1 )
    t(1) = dtd2*bx0;
    t(2) = dtd2*pb(2);
    t(3) = dtd2*pb(3);
    opts = 1 + t(1)^2 + t(2)^2 + t(3)^2;
    s(1) = 2*t(1)/opts;
    s(2) = 2*t(2)/opts;
    s(3) = 2*t(3)/opts;
  end

% Now step particle positions and velocities
  
   if( ipc==1 )
     px(m,1) = 0.5*( px(m,1) + px(m,2) ) + dt*pv(m,1,2);
   else
     px(m,2) = px(m,1) + dt*pv(m,1,2);
   end

   if( iboris==1 )
     if( ipc==1 )
       pvm(1) = 0.5*( pv(m,1,1) + pv(m,1,2) ) + dtd2*pe(1);
       pvm(2) = 0.5*( pv(m,2,1) + pv(m,2,2) ) + dtd2*pe(2);
       pvm(3) = 0.5*( pv(m,3,1) + pv(m,3,2) ) + dtd2*pe(3);
     else
       pvm(1) = pv(m,1,1) + dtd2*pe(1);
       pvm(2) = pv(m,2,1) + dtd2*pe(2);
       pvm(3) = pv(m,3,1) + dtd2*pe(3);
     end

     pvs(1) = pvm(1) + pvm(2)*t(3) - pvm(3)*t(2);
     pvs(2) = pvm(2) + pvm(3)*t(1) - pvm(1)*t(3);
     pvs(3) = pvm(3) + pvm(1)*t(2) - pvm(2)*t(1);

     pvp(1) = pvm(1) + pvs(2)*s(3) - pvs(3)*s(2);
     pvp(2) = pvm(2) + pvs(3)*s(1) - pvs(1)*s(3);
     pvp(3) = pvm(3) + pvs(1)*s(2) - pvs(2)*s(1);

     if( ipc==1 )
       pv(m,1,1) = pvp(1) + dtd2*pe(1);
       pv(m,2,1) = pvp(2) + dtd2*pe(2);
       pv(m,3,1) = pvp(3) + dtd2*pe(3);
     else
       pv(m,1,2) = pvp(1) + dtd2*pe(1);
       pv(m,2,2) = pvp(2) + dtd2*pe(2);
       pv(m,3,2) = pvp(3) + dtd2*pe(3);
     end
   else
     if( ipc==1 )
       pv(m,1,1) = 0.5*( pv(m,1,1) + pv(m,1,2) ) + ...
                  dt*( pe(1) + pv(m,2,2)*pb(3) - pv(m,3,2)*pb(2) );
       pv(m,2,1) = 0.5*( pv(m,2,1) + pv(m,2,2) ) + ...
                  dt*( pe(2) + pv(m,3,2)*bx0   - pv(m,1,2)*pb(3) );
       pv(m,3,1) = 0.5*( pv(m,3,1) + pv(m,3,2) ) + ...
                  dt*( pe(3) + pv(m,1,2)*pb(2) - pv(m,2,2)*bx0   );
     else
       pv(m,1,2) = pv(m,1,1) + dt*( pe(1) + pv(m,2,2)*pb(3) - pv(m,3,2)*pb(2) );
       pv(m,2,2) = pv(m,2,1) + dt*( pe(2) + pv(m,3,2)*bx0   - pv(m,1,2)*pb(3) );
       pv(m,3,2) = pv(m,3,1) + dt*( pe(3) + pv(m,1,2)*pb(2) - pv(m,2,2)*bx0   );
     end
   end
  
end

% Fix points that have gone beyond the boundaries

for m= 1:nm
  if( px(m,ipc)>=lsim )
    px(m,ipc) = px(m,ipc) - lsim;
  elseif( px(m,ipc)<0 )
    px(m,ipc) = px(m,ipc) + lsim;
  end
end
