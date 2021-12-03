

zimboras_gr = function()
{
  a = gfc("L4+C3+L4")
  a$target = 7
  a = join_link(a,gfc("L4"))

  return(a)
}

correct_cut_square_unit = function(N = 4)
{
  out = ringNcut(N,c(1,1 + N%/%2), HANDLES = F)
  out$re_coord[[2]] = c(1,N)

  return(out)
}

eigen_uncertainty = function(gr , sample = 20)
{
  eigen_list = lapply(  rep(list(gr$mat), sample), eigen)

  value_list = sapply( eigen_list, "[[", "values")

  mean  = apply(value_list,1,mean)
  unc   = apply(value_list,1,sd)

  print(value_list[1,])

  print( matrix(c(mean,unc), ncol = 2))


}
