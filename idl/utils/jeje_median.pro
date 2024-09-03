function jeje_median,a

  ind = sort(a)
  med = a(ind(n_elements(a)/2))
  return,med
  
end
