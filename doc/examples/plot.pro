;listing00
pro test 
  ; specifying parameters
  gstdev = 1.5
  outfile = 'idltest.nc'

  ; executing the simulation
  spawn, [ $
    'python', 'parcel.py', '--outfile', outfile, '--gstdev', string(gstdev) $
  ], /noshell

  ; opening the output file
  nc = ncdf_open(outfile, /nowrite)

  ; plotting initial and final wet spectra
  ncdf_diminq, nc, ncdf_dimid(nc, 'radii'), ignore, nr
  ncdf_diminq, nc, ncdf_dimid(nc, 't'    ), ignore, nt

  !P.MULTI=[0,2,1]
  ncdf_varget, nc, 'radii_dr', radii_dr
  ncdf_varget, nc, 'radii_rl', radii_rl
 
  foreach it, [0, nt-1] do begin
    ncdf_varget, nc, 'radii_m0', radii_m0, count=[nr, 1], offset=[0, it]
    plot,                                $
      xtitle='particle wet radius [um]', $
      ytitle='[mg-1 um-1]',              $
      psym=10, /xlog, /ylog,             $
      (radii_rl+radii_dr/2)*1e6, (radii_m0/1e3)/(radii_dr*1e6)
  endforeach
end
;listing01
