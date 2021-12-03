# linear model extrapolation

# It's up to the user to set the right global parameter to do the linear extrapolation
# Those methods pick a set of ~15 contiguous size graph and tries a lm on their behavour

# Best suited for GLOBAL_FIRST_MAX set to 1

line_t_lm = function(mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  range = 5:20

  gr_vec = lapply(range,line)
  #no need to optimize

  t_vec = sapply(gr_vec, t_p_max_func)


  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )
}

ring_t_lm = function(mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  #odd ring graph tend to be opt at first maximum
  range =5:20

  gr_vec = lapply(range,ring)
  gr_vec = lapply(gr_vec,optimize_gr)

  t_vec = sapply(gr_vec, t_p_max_func)

  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )
}

ringNh_t_lm = function(mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  range =5:20

  gr_vec = lapply(range,ringNhandle)
  gr_vec = lapply(gr_vec,optimize_gr)

  t_vec = sapply(gr_vec, t_p_max_func)


  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )
}

ringNc_t_lm = function(mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  range =5:20

  gr_vec = lapply(range, ringNcut , HANDLES = F)
  gr_vec = lapply(gr_vec,optimize_gr)

  t_vec = sapply(gr_vec, t_p_max_func)


  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )
}

#get linear model of time of arrival of first maxima ring couples, with first ring fixed and the secon widening
ringCouple_t_lm = function(handles = 0,first_ring = 4, op = "/", mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  range =5:30
  base_str = "La+CbOCN+La"
  base_str = gsub("a", handles,base_str)
  base_str = gsub("b", first_ring,base_str)
  base_str = gsub("O", op,base_str)

  str_vec = sapply(range,nsub,codes = base_str)

  gr_vec = lapply(str_vec,gfc)
  gr_vec = lapply(gr_vec,optimize_gr)

  t_vec = sapply(gr_vec, t_p_max_func)


  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )

}

#get linear model of time of arrival of first maxima for handled ring, withring widening
handles_t_lm = function(handles = 2, mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  range =5:30
  base_str = "La+CN+La"
  base_str = gsub("a", handles,base_str)

  str_vec = sapply(range,nsub,codes = base_str)

  gr_vec = lapply(str_vec,gfc)
  gr_vec = lapply(gr_vec,optimize_gr)

  t_vec = sapply(gr_vec, t_p_max_func)


  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )
}

#get linear model of time of arrival of first maxima for handled ring, with handles widening
cycle_t_lm = function(cycle = 2, mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  range =0:20
  base_str = "La+CN+La"
  base_str = gsub("N", cycle,base_str)

  str_vec = sapply(range,gsub, pattern = "a",x = base_str)

  gr_vec = lapply(str_vec,gfc)
  gr_vec = lapply(gr_vec,optimize_gr)

  t_vec = sapply(gr_vec, t_p_max_func)


  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )
}

#get linear model of time of arrival of first maxima for chains composed by unit
chain_t_lm = function(unit, space = 0, mode = "dist", TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  #trend stabilizes for larger N
  range =1:20

  gr_vec = lapply(range,chain, gr = unit, space = space)
  gr_vec = lapply(gr_vec,optimize_gr, mode = "diag")

  t_vec = sapply(gr_vec, t_p_max_func)

  #print(t_vec)


  x = get_list_x(gr_vec,mode)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ x) )

}


#meh
ringNh_p_bigN_lm = function(  TC = 1 , FIRST = T)
{
  pushTC(TC)
  pushFC(FIRST)

  range = 20:35

  gr_vec = lapply(range,ringNhandle)
  gr_vec = lapply(gr_vec,optimize_gr)

  t_vec = sapply(gr_vec, p_max_func)

  #pop old value
  popTC()
  popFC()

  return( lm( t_vec ~ range) )
}

list_t_lm = function(gr_vec, mode = "dist")
{
  t_vec = sapply(gr_vec, t_p_max_func)

  x = get_list_x(gr_vec,mode)

  return( lm( t_vec ~ range) )
}

list_p_lm = function(gr_vec, mode = "dist")
{
  p_vec = sapply(gr_vec, p_max_func)

  x = get_list_x(gr_vec,mode)

  return( lm( p_vec ~ range) )
}

plot_handles_lm_coeff = function(range = 0:10)
{
  lm_list = lapply(range,handles_t_lm)

  plot_lm_list(lm_list, range = range, x_lab = "handles size")
}

plot_cycle_lm_coeff = function(range = 3:13)
{
  lm_list = lapply(range,cycle_t_lm)

  plot_lm_list(lm_list, range = range, x_lab = "central ring size")
}

plot_couple_lm_coeff = function(range = 3:10)
{
  lm_list = lapply(range,ringCouple_t_lm, handles = 0)

  plot_lm_list(lm_list, range = range, x_lab = "first cycle size")
}

plot_chain_lm_coeff = function(unit, range = 0:10)
{
  lm_list = lapply(range,chain_t_lm, unit = unit)

  plot_lm_list(lm_list, range = range, x_lab = "space btw units")
}

plot_lm_list = function( lm_list, range, x_lab)
{
    inter_vec = sapply(lm_list, "[[", "coefficients")[1,]
  slope_vec = sapply(lm_list, "[[", "coefficients")[2,]
  line_coeff = line_t_lm()$coefficients

  plot(-10,-10,
       ylim = c(-1,3), xlim= range( c(-1,range) ),
       xlab = x_lab, ylab="coeff")
  title(main ="")

  legend("topright",
          legend = c("intercept", "slope"),
          col = c(2,3,1), lty = c(1,1,-1),pch = c(1,1,19))

  lines(range, inter_vec, col = 2)
  points(range, inter_vec, col = 2)
  points(-1,line_coeff[1], col = 2, pch = 19)

  lines(range, slope_vec, col = 3)
  points(range, slope_vec, col = 3)
  points(-1,line_coeff[2], col = 3, pch = 19)

  grid()
}
