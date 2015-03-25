function b = step_b( b, e, ni, nim, dtddx2, ipc )

% = 1/2 for predictor/corrector step

if( ipc==1 )
  for i= 2:nim
    b(i,2,1)= 0.5*( b(i,2,1) + b(i,2,2) ) + ( e(i+1,3) - e(i-1,3) )*dtddx2;
    b(i,3,1)= 0.5*( b(i,3,1) + b(i,3,2) ) - ( e(i+1,2) - e(i-1,2) )*dtddx2;
  end
else
  for i= 2:nim
    b(i,2,2)= b(i,2,1) + ( e(i+1,3) - e(i-1,3) )*dtddx2;
    b(i,3,2)= b(i,3,1) - ( e(i+1,2) - e(i-1,2) )*dtddx2;
  end
end

b( 1,2:3,ipc)= b(nim,2:3,ipc);
b(ni,2:3,ipc)= b(  2,2:3,ipc);


end
