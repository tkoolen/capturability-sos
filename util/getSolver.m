function solver = getSolver()

ok_mosek = checkDependency('mosek');
ok_sedumi = checkDependency('sedumi');

if ~ok_sedumi && ~ok_mosek
  error('You need either MOSEK or SeDuMi installed to use this function.');
end

% Choose sdp solver
if ok_mosek
  solver = @spot_mosek;
else
  solver = @spot_sedumi;
end

end