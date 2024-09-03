pro read_acc_story,outputdir,type,acc

  if type eq 'acc' then $
     filename = outputdir + '/acc_story.dat' $
  else $
     filename = outputdir + '/wind_story.dat'
  
  openr,11,filename
  r=(n=0L)
  readu,11,r,n,r
  time = dblarr(n)
  readu,11,r,time,r
  ; count histories
  ng = 0L
  while not eof(11) do begin 
     ng++
     readu,11,r,time,r
  endwhile
  close,11
  
  acc = {ng:0L,ntime:0L,time:dblarr(n),acc:dblarr(ng,n)}

  openr,11,filename
  r=(n=0L)
  readu,11,r,n,r
  acc.ntime = n
  acc.ng = ng
  time = dblarr(n)
  readu,11,r,time,r
  acc.time = time
  ; read histories
  a = dblarr(n)
  for ig = 0L,ng-1L do begin 
     readu,11,r,a,r
     acc.acc(ig,*) = a
  endfor
  close,11

end
