



pro mainfun6

  ; define k values before running the program
  ; k represents for band interval; when k=0, no band interval, bands are collected continuously; k=1 meaning b1, b3, b5...
  ; k< (l-n)/n, l is number of bands, n is number of input spectra
  ; example: k=[0,3,7,15,31,63,127,255], usually k is defined as the first small numbers
  k=[0,1,3,7]
  
  ; input is a spectral library file in ENVI
  envi_select, title='Input endmember spectral library', fid=sfid, pos= pos, dims=sdims
  if sfid eq -1 then return
  envi_file_query, sfid, ns=sns, nl=snl, file_type= sfile_type, data_type=sdata_type, wl=swl, spec_names=sspec_names, FNAME=sCoreFileName

  if (sfile_type ne 4) then begin
    ret=dialog_message('None Spectral Library format')
    return
  endif

  wavelength=swl
  data=ENVI_get_data(fid=sfid, dims=sdims, pos=pos)


  E=transpose(data)
  size_E=size(E,/dimensions)
  l=size_E(0)
  r=size_E(1)


  k0=k(0)
  k1=n_elements(k)

  angle=dblarr(k1,r)
  E1=dblarr(l,l)
  ;  ; normalize the E matrix
  ;  for i=0,l-1 do begin
  ;    E[i,* ]=E[i, *]/norm(E[i,*],/double,lnorm=2)
  ;  endfor

  ;
  ;  for i=0,r-l do begin
  ;    angle(0,i+(l-1)/2)=nsolidangle(E(*,i:(i+l-1)))
  ;  endfor

  for m=0,k1-1 do begin
    k0=k(m)
    for j=0,r-(l-1)*(k0+1)-1 do begin

      for i=0,l-1 do begin
        E1(*,i)=E(*,j+i*(k0+1))
      endfor
      angle(m,j+((l-1)*(k0+1))/2)=nsolidangle(E1)
      
    endfor

  endfor
  
  rankings=dblarr(k1,r)
  band_index=dblarr(k1,r)
  for i=0,k1-1 do begin

    rankings(i,*)=wavelength[reverse(sort(angle(i,*)))]
    band_index(i,*)=reverse(sort(angle(i,*)))
    ;  bands_output=transpose(rankings[0:100])
    ;  bands_output=transpose(rankings)
  endfor
  
  print, angle
  print,'band ranking wavelength='
  print,rankings
  print,'band index='
  print, band_index+1
  
  
  ; output the result in text file
  out_file1=dialog_pickfile(title='Output band ranking File',default_extension='txt',file='band_ranking')
  openw,lun1,out_file1,/get_lun
  printf,lun1,rankings, FORMAT='('+strtrim(k1,2)+'E)' 
  ;  for j=0,k1-1 do printf, lun, rankings(j,*)
  free_lun,lun1
  ;
  
  out_file2=dialog_pickfile(title='Output band index File',default_extension='txt',file='band_index')
  openw,lun2,out_file2,/get_lun
  printf,lun2,band_index+1, FORMAT='('+strtrim(k1,2)+'I)'
  ;  for j=0,k1-1 do printf, lun, rankings(j,*)
  free_lun,lun2
  
  out_file3=dialog_pickfile(title='Output angle File',default_extension='txt',file='angle')
  openw,lun3,out_file3,/get_lun
  printf,lun3,angle, FORMAT='('+strtrim(k1,2)+'G)' 
  close,/all




end


function f_angle, E,g

  size_E=SIZE(E,/dimensions)
  l=size_E(0)


  for i=0,l-1 do begin
    E[i,* ]=E[i, *]/norm(E[i,*],/double,lnorm=2)
  endfor
  ;--------------------------------------------------------------------------------------------------------------------
  ; this is the beginning
  A=dblarr(2,l)
  for i=0,l-2 do begin
    A[*,i]=[sin(g(i)),cos(g(i))]
  endfor



  V=fltarr(1,l)

  t=1
  for i=0,l-2 do begin
    t=t*A(0,i)
  endfor
  V(l-1)=t

  for i=1,l-2 do begin
    t=1
    for j=1,i do begin
      t=t*A(0,j-1)
    endfor
    V(i)=t*A(1,i)
  endfor

  V(0)=A(1,0)




  jacobia=1
  for i=0,l-3 do begin
    jacobia=jacobia*((A(0,i))^(l-i-2))
  endfor




  t=0

  for i=0,l-2  do begin
    for j=i+1,l-1 do begin


      xishu=double(total(E(i,*)*E(j,*)))
      t=(t+ xishu* V(i)*V(j))
    endfor
  endfor

  ;  print,  xishu


  ft=(1+2*t)^(l/2)

  ;  jj=jacobian(g(0),g(1),g(2),g(3))

  f=jacobia/ft



  return, f

end


function nsolidangle,E

  N=1000


  size_E=SIZE(E,/dimensions)
  l=size_E(0)

  ;
  ;  for i=0,l-1 do begin
  ;    E[i,* ]=E[i, *]/norm(E[i,*],/double,lnorm=2)
  ;  endfor
  ;
  ; using Monte Carlo to calculate the integral

  ; g is created for 'ndimensionangle.pro' function

  g=randomn(seed,l-1)

  Q=(dindgen(1,l-1)+1)*0.1
  R=10*!pi*tan(Q)


  N1=indgen(N,1)+1
  N2=intarr(N,l-1)
  for i=0,l-2 do begin
    N2(*,i)=N1
  endfor

  R1=dindgen(N,l-1)
  for i=0,N-1 do begin
    R1(i,*)=R
  endfor

  Ran=transpose(N2*R1)


  Ran=Ran mod (!pi/2)


  p=dblarr(1,N)

  ; call ndimensionanlge.pro function by N times
  for i=0,N-1 do begin

    for j=0, l-2 do begin
      g(j)=Ran(j,i)

    endfor

    p(i)=f_angle(E,g)

  endfor
  ; get the angle
  pp=total(p)
  angle=abs(determ(E))*((!pi/2)^(l-1))*(double(1)/double(N))*pp
  ;    angle=double(abs(determ(E))*((!pi/2)^(l-1))*(1/N)*(pp))

  return, angle

end


