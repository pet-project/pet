clear all;
close all;

%
% Input parameters
%

% nmd8= the number of particles divided by 8.

nmd8 = 16;

% nid= the number of data grid points.  These are numbered from 
%  2 to nid+1.  The dimension ni=nid+2 holds 2 extra "buffer" points 
%  numbered 1 and ni=nid+2.  These two points wrap around periodically
%  so that rho(1)=rho(ni-1) and rho(ni)=rho(2).

nid = 16;

% thetakb = angle (in degrees) between B_0 and the x axis

thetakb= 0;

% lsim = length of simulation in units of c/omega_pi

lsim= 2;

% nsamplev= the number of sample velocities used for creating a 
%  Maxwellian velocity distribution function for v_par.

nsamplev = 500;

% beta defined using parallel temperature 

beta_par = 0;

% trat = T_perp / T_par - ratio of perpendicular to parallel temperatures

trat= 2;
      
% fmax.  vth*fmax is the largest velocity represented by the sample
%  velocity array samplev.  Since the Maxwellian velocity distribution
%  falls off fairly rapidly, it is not necessary for this to be very 
%  large.  Somewhere between 3 and 4 is fine.

fmax = 3;

% i_quiet_start = 1 for quiet start

i_quiet_start = 0;

% vzpert = velocity perturbation in y direction

vzpert= 1.e-3;

% vlpert = amplitude of left hand linear wave

vlpert = 0.e-3;

% vrpert = amplitude of right hand linear wave

vrpert = 0.e-3;

% iboris = 1 to use Boris mover

iboris= 1;

% dt= time step (normalized to cyclotron time).

dt = .01;

% tsim = time of simulation.

tsim = 30;  

% tfldout= number of time steps between field outputs

tfldout = .05;

% tftout = time between outputs to array for fourier transforms

tftout = .02;

% ipartinc= increment between plotted particles
%   = 1 to plot all points
%   = larger number to speed up code

ipartinc = 1;

% tpause = time to pause (in s) between field outputs

tpause = 0.2;

% i_subtract_thermal_energy = 1 to subtract off thermal part of kinetic
% energy

i_subtract_thermal_energy = 1;

%
% Subsidiary parameters (depending on input paramters above)
%

nm = nmd8*8;               % total number of particles
nmd2 = nm/2;
ni = nid+2;                % total number of grid points including the buffer zone
pweight= nid/nm;           % particle weight
nim = ni - 1;
dx = lsim/nid;             % distance between grid points
dx2 = 2*dx;
dtd2 = dt/2;
dtddx2 = dt/dx2;
k1 = 2*pi/lsim;            % k of longest wavelength mode
nt = tsim/dt;              % number of time steps in simulation
nfldout = tfldout/dt;      % number of time steps between plots
nftout = tftout/dt;        % number of time steps between Fourier transform outputs
iift = 2 + nid/4;          % grid point for Fourier transform data
t_par = beta_par/2;        % parallel temperature
t_perp = t_par * trat;     % perp temperature
vth_par = sqrt(t_par);                      % parallel thermal velocity
vth_perp = sqrt(t_perp);                    % perp thermal velocity
fmax_perp = sqrt(2)*fmax;                   % fmax for v_perp
vmaxy_plot = max([3*vth_perp ; vzpert ; vlpert ; vrpert]);    % maximum v for v space plots
expperp = exp(-fmax_perp^2/2);              % useful quantity for Maxwellian distribution
bx0= cos(thetakb*pi/180);                   % x component of backgroup magnetic field
by0= sin(thetakb*pi/180);                   % y  component of backgroup magnetic field


px = zeros(nm,2);          % particle position, 2 times
pv = zeros(nm,3,2);        % particle velocity, 3 components, 2 times

x = zeros(ni,1);           % position
b = zeros(ni,3,2);         % B, 3 components, 2 times
e = zeros(ni,3);           % E, 3 components
cur = zeros(ni,3);         % J, 3 components
curp = zeros(ni,3);        % J_p (proton current), 3 components
ve = zeros(ni,3);          % v_e (electron velocity), 3 components
ne = zeros(ni,1);            % n_e

% samplev= array of sample velocities used for creating a Maxwellian
%  velocity distribution function.

samplev = zeros(nsamplev,1);

% Define x positions. x=0 and 1 at grid points 1.5 and ni-.5.
% This array is helpful for plotting rho 
%  and e.  Note that x=0 inbetween the first (buffer)
%  grid point and the first (i=2) data point.

for i= 1:ni
  x(i)= (i-1.5)*dx;
  b(i,1,1)= bx0;           % array not really needed since b_x = constant
  b(i,2,1)= by0;
  b(i,3,1)= 0;
end
cur(:,1)= 0;               % cur_x = 0 in 1D

% Set up array samplev with probabilities associated with certain 
%  velocities.

[dv,samplev] = setmaxw(fmax,nsamplev);

% Define initial particle positions (uniform distribution in x with
%  sinusoidal perturbation) and velocities.  Here pi2=2*pi in radians.

pi2 = 2*pi;
if( i_quiet_start==1 )
  pid2 = pi/2;
  for md8= 1:nmd8
    x_value = (md8-.5)*lsim/nmd8;
    r_par = nrev(md8,2);
    vthpart_par = vth_par * ...
                    getmaxw(r_par,samplev,nsamplev,dv);
    r_perp = nrev(md8,3);
    vthpart_perp = vth_perp * ...
            sqrt( -2*log( 1 - r_perp*( 1 - expperp ) ) );
    phase0 = 2*pi*nrev(md8,5);
    vzpertpart= vzpert*sin(k1*x_value);
    m1m = (md8-1)*8;
    for m = m1m+1:m1m+8
      px(m,1) = x_value;
    end
    for m = m1m+1:m1m+4
      pv(m  ,1,1) = -vthpart_par;
      pv(m+4,1,1) = +vthpart_par;
    end
    pv(m1m+1,2,1) = vthpart_perp*cos(phase0);
    pv(m1m+1,3,1) = vthpart_perp*sin(phase0);
    pv(m1m+2,2,1) = vthpart_perp*cos(phase0+pid2);
    pv(m1m+2,3,1) = vthpart_perp*sin(phase0+pid2);
    pv(m1m+3,2,1) = vthpart_perp*cos(phase0+2*pid2);
    pv(m1m+3,3,1) = vthpart_perp*sin(phase0+2*pid2);
    pv(m1m+4,2,1) = vthpart_perp*cos(phase0+3*pid2);
    pv(m1m+4,3,1) = vthpart_perp*sin(phase0+3*pid2);
    for m = m1m+1:m1m+4
      pv(m+4,2,1) = pv(m,2,1);
      pv(m+4,3,1) = pv(m,3,1);
    end
    for m = m1m+1:m1m+8
      pv(m,3,1)= pv(m,3,1) + vzpertpart;
    end
  end
else
  for m = 1:nmd2
    px(m,1) = (m-.5)*lsim/nmd2;
    vzpertpart= vzpert*sin(k1*px(m,1));
    r_par = rand;
    pv(m,1,1) = vth_par * ...
                      getmaxw(r_par,samplev,nsamplev,dv);
    r_perp = rand;
    vperp = vth_perp * ...
          sqrt( -2*log( 1 - r_perp*( 1 - expperp ) ) );
    theta = pi2*rand;
    pv(m,2,1) = vperp * cos(theta);
    pv(m,3,1) = vperp * sin(theta) + vzpertpart;

    px(m+nmd2,1)   =  px(m,1);
    pv(m+nmd2,1,1) = -pv(m,1,1);
    pv(m+nmd2,2,1) =  pv(m,2,1);
    pv(m+nmd2,3,1) =  pv(m,3,1);
  end
end 

%Fix up this section:

% Add in left and right hand perturbations

% for ilr = -1:2:1      % - sign for left hand mode
%   ilr
%   if( ilr==-1 )
%     vlrpert= vlpert;
%   else
%     vlrpert= vrpert;
%   end
%   k1s = k1^2;
%   k14 = k1s^2;
%   omegas = ...
%   omegalr = sqrt( omegas );
%   vpy0 = vlrpert;
%   vpz0 = ...
%   dby0 = ...
%   dbz0 = ...
%   for m = 1:nm
%     pv(m,2,1) = pv(m,2,1) + ...
%     pv(m,3,1) = pv(m,3,1) + ...
%   end
%   for i = 1:ni
%     b(i,2,1) = b(i,2,1) + ...
%     b(i,3,1) = b(i,3,1) + ...
%   end  
% end

% Now that the perturbations are added, initialize second copy of particle and b arrays
b(:,:,2)= b(:,:,1);
px(:,2)= px(:,1);
pv(:,:,2)= pv(:,:,1);

% Evaluate subsidiary quantities
[ne, curp] = from_particles(ne,curp,ni,x(1),dx,px,pv,nm,pweight);
[cur,ve,e] = from_fields(cur,ve,e,b,bx0,curp,ne,ni,nim,dx2);

make_two_figures;

it = 1;
t = 0
tarray(1) = 0;
pen1 = 0;
[ben(it),pen(it),ten(it)] = energy(it,b,ni,nid,by0,pv,nm,pweight,pen1,i_subtract_thermal_energy );
pen1 = pen(1);
if( i_subtract_thermal_energy==1 ) 
  pen(1) = 0;
  ten(1) = ben(1);
end
ift = 1;
tft(1)= 0;
bya(1)= b(iift,2,2);
bza(1)= b(iift,3,2);

% Do simulation for nt-1 time steps (nt times). Set initial time t to 0.

b2max= -1;
b3max= -1;
psb2max= -1;
bamax= -1;
for it= 2:nt+1

    for ipc= 1:2           % ipc = 1/2 for predictor/corrector

        [px,pv] = step_particles_boris( ...
            px, pv, nm, e, b, ni, nim, lsim, x(1), dx, bx0, dt, dtd2, iboris, ipc );
        b = step_b( b, e, ni, nim, dtddx2, ipc );
     
        if( ipc == 1 )
            if( it==2 )
                px(:,    1) = 0.5 .* ( px(:,    1) + px(:,    2) );
                pv(:,  :,1) = 0.5 .* ( pv(:,  :,1) + pv(:,  :,2) );
                b (:,2:3,1) = 0.5 .* (  b(:,2:3,1) + b (:,2:3,2) );
            end  
            [px(:,    1),px(:,    2)] = switch_aa_aa(px(:,    1),px(:,    2) );
            [pv(:,  :,1),pv(:,  :,2)] = switch_aa_aa(pv(:,  :,1),pv(:,  :,2) );
            [ b(:,2:3,1), b(:,2:3,2)] = switch_aa_aa( b(:,2:3,1), b(:,2:3,2) );
        end
     
        [ne, curp] = from_particles(ne,curp,ni,x(1),dx,px,pv,nm,pweight);
        [cur,ve,e] = from_fields(cur,ve,e,b,bx0,curp,ne,ni,nim,dx2);
    end

    [ben(it),pen(it),ten(it)] = energy(it,b,ni,nid,by0,pv,nm,pweight,pen1,i_subtract_thermal_energy );
    t= t + dt;
    tarray(it) = t;
  
    if( mod(it-1,nftout)==0 ) 
        ift= ift + 1;
        tft(ift) = t;
        bya(ift) = b(iift,2,2);
        bza(ift) = b(iift,3,2);
        bamax = max([bamax ; abs(bya(ift)) ; abs(bza(ift))]);
    end

  % Make plots.

  if( mod(it-1,nfldout)==0 )

    t
    
    figure(1)
 
    subplot(4,2,1)
    plot(x,b(:,2,2),'b')
    xlim([x(1) x(ni)])
    b2max = max([ b2max ; abs(b(:,2,2)) ]);
    ylim([-b2max b2max])
    xlabel('x')
    ylabel('B_y')
    title(strcat('time=',num2str(t)));

    subplot(4,2,3)
    plot(x,b(:,3,2),'r')
    xlim([x(1) x(ni)])
    b3max = max([ b3max ; abs(b(:,3,2)) ]);
    ylim([-b3max b3max])
    xlabel('x')
    ylabel('B_z')
    
    v_array = -vmaxy_plot + ((1:21)-1)*2*vmaxy_plot/20;
    subplot(4,2,2)
    hist(pv(:,1,2),v_array)
    xlabel('v_{//}')
    ylabel('N_v')

     
    subplot(4,2,4)
    hist(pv(:,2,2),v_array)
    xlabel('v_y')
    ylabel('N_v')
    
    subplot(4,1,3)
    kby = ifft(b(2:nim,2,2));
    pskby = kby.*conj(kby);
    pskby(2:nid/2) = 2*pskby(2:nid/2);
    kbz = ifft(b(2:nim,3,2));
    pskbz = kbz.*conj(kbz);
    pskbz(2:nid/2) = 2*pskbz(2:nid/2);
    psb2max = max([ psb2max ; pskby(1:nid/2+1) ; pskbz(1:nid/2+1) ]);
    plot((0:nid/2),pskby(1:nid/2+1),'b',(0:nid/2),pskbz(1:nid/2+1),'r')
    xlim([0 nid/2])
    ylim([0 psb2max])
    xlabel('Mode Number')
    ylabel('Power')

    if( ift>=8 )
      subplot(4,1,4)   
      npow2 = 2^(floor(log2(ift)));
      npow2d2 = npow2/2;
      ba = bya(1:npow2) + 1i*bza(1:npow2);
      fba = ifft(ba);
      sfba = fftshift(fba);
      psfba = sfba.*conj(sfba);
      [psfbamax,imax] = max(psfba);
      i1 = find(psfba>0.05*psfbamax,1,'first');
      i1 = max([i1-1 1]);
      i2 = find(psfba>0.05*psfbamax,1,'last' );
      i2 = min([i2+1 npow2]);
      tlen = npow2*tftout;
      omega0 = 2*pi*(1/tlen);
      omega = -(npow2d2-1)*omega0 + ((1:npow2)-1)*omega0;
      plot(omega,psfba)
      xlim([omega(i1) omega(i2)])
      xlabel('\omega')
      ylabel('Power')
    end
    
    figure(2)
    
    subplot(3,3,1)
    plot(bya(ift),bza(ift),'.')
    bamaxp= bamax*1.2;
    xlim([-bamaxp bamaxp])
    ylim([-bamaxp bamaxp])
    xlabel('bya')
    ylabel('bza')

    subplot(3,3,2)
    plot(pv(1:ipartinc:nm,1,2),pv(1:ipartinc:nm,3,2),'.')
    xlim([-vmaxy_plot vmaxy_plot])
    ylim([-vmaxy_plot vmaxy_plot])
    xlabel('v_x')
    ylabel('v_z')
    title(strcat('time=',num2str(t)));
    
    subplot(3,3,3)
    plot(pv(1:ipartinc:nm,2,2),pv(1:ipartinc:nm,3,2),'.')
    xlim([-vmaxy_plot vmaxy_plot])
    ylim([-vmaxy_plot vmaxy_plot])
    xlabel('v_y')
    ylabel('v_z')
    
    subplot(6,1,3)
    plot(tft,bya)
    xlim([0 tsim])
    ylim([-bamax bamax])
    xlabel('t')
    ylabel('bya')
    
    subplot(6,1,4)
    plot(tft,bza)
    xlim([0 tsim])
    ylim([-bamax bamax])
    xlabel('t')
    ylabel('bza')
       
    subplot(6,1,5)
    semilogy(tarray(2:it),ben(2:it))
    xlim([0 tsim])
    xlabel('t')
    ylabel('E_B')
    
    subplot(6,1,6)
    plot(tarray,pen,'r',...
      tarray,ben,'b',...
      tarray,ten,'k')
    xlim([0 tsim])
    xlabel('t')
    ylabel('E_B, E_K, E_{tot}')
    
    pause(tpause)
  end;

end;  % end of 'for it='
