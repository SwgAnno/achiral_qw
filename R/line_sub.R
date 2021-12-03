
# plot how best performance varies as the central links of a line are divided,
# forming an increasing size central ring

plot_line_to_ring = function( L)
{

  plot(-10,-10,ylim = c(0,1.3), xlim= c(0,L),xlab = "subbed links", ylab=expression(p["max"]))
	legend("topleft", title = paste("L",L," subs",sep=""), ncol = 2,
         legend = c("odd sub", "even_sub"," '' nophi", " '' nophi"),
         col = c(2:3,2:3),lty= c(1,1,3,3), pch = c(19,19,-1,-1))
  grid()

  gr_list = list()

 #always even

  base_str = "La/Cb/Lc"

  for(i in 1:(L-1) )
  {
    dx = (L-1-i)%/%2 +1
    sx = dx+ (L-1-i)%%2

    #print(c(sx,dx,(L-i)))

      gr_str = gsub("a",sx,base_str)
      gr_str = gsub("c",dx,gr_str)
      gr_str = gsub("b",2*i,gr_str)

      gr_list[[2*i -1]] = gfc(gr_str)

      gr_str = gsub("a",sx,base_str)
      gr_str = gsub("c",dx,gr_str)
      gr_str = gsub("b",2*i+1,gr_str)
      gr_list[[2*i]] = gfc(gr_str)
  }


  gr_list_odd = gr_list[2*(1:(L-1))]
  gr_list_even = gr_list[2*(1:(L-1)) -1]

  opt_and_plot_list(gr_list_odd,mode = "range", range =1:(L-1),nophi = T, pch = 19, col = 2)
  opt_and_plot_list(gr_list_even,mode = "range", range =1:(L-1),nophi = T, pch = 19, col = 3)

  return(sapply(gr_list,"[[","code"))

}

#plot how best performance varies as a base unit gets added increasing length handles
plot_add_handles = function(gr, range =1:10)
{
  plot(-10,-10,ylim = c(0,1.3), xlim= c(0, max(range)),xlab = "handle size", ylab=expression(p[max]))
	legend("topleft",ncol = 2,
         legend = c("both", "base","dx/sx" ,"1sx +dx"),
         col = c(2,3,5,6), pch = 19)
  grid()

  points( 0, p_max_func(optimize_gr(gr)), pch = 19, col = 3)

  #####
  gr_list = lapply(range, add_handles, base = gr, mode = "both")

  print(sapply(gr_list,"[[", "code"))

  opt_and_plot_list(gr_list,mode = "range", range =range,nophi = TRUE, col = 2, pch = 19)
  print("pappa")


  #####
  gr_list = lapply(range, add_handles, base = gr, mode = "dx")

  opt_and_plot_list(gr_list,mode = "range", range =range,nophi = TRUE, col = 4, pch = 19)

  #####
  gr_list = lapply(range, add_handles, base = gr, mode = "sx")

  opt_and_plot_list(gr_list,mode = "range", range =range,nophi = TRUE, col = 5, pch = 19)

  #####
  gr_list = lapply(range-1, add_handles, base = gr, mode = "fixr", fix = 1)

  opt_and_plot_list(gr_list,mode = "range", range =range-1,nophi = TRUE, col = 6, pch = 19)



}

#plot how best performance varies as a substitution unit gets inserted inside a line
plot_slide_into_line = function(L, in_gr)
{
  gr_list = list()

  tot_dist = L-2 + distance(in_gr)

  for(i in 1:(L-1))
  {
    sx = i
    dx = L-i

    #add left&right appendix
    temp_gr = join_nolink(line(sx),in_gr)
    temp_gr = join_nolink(temp_gr,line(dx))

    gr_list[[i]] = temp_gr

  }

  max_vec_nophi = sapply(gr_list,p_max)

  gr_list = lapply(gr_list,optimize_gr)

  max_vec = sapply(gr_list,p_max)

  plot(-10,-10,ylim = c(0,1.3), xlim= c(1,L-1),xlab = "subbed link", ylab="p max")
	legend("topleft", title = paste("L",L," with sliding gr",sep=""),
         legend = c("gr", " '' nophi", paste("L",tot_dist,sep="")),
         col = c(2,2,3),lty= c(1,3,3), pch = c(19,-1))
  grid()

  sites = 1:(L-1)

  lines( sites, max_vec, pch = 19, col = 2)
  lines( sites, max_vec_nophi, pch = 19, col = 2, lty = "dotted")
	points( sites, max_vec, pch = 19, col = 2)

	abline( h = p_max(line(tot_dist)),col=3, lty = "dotted")

  return(sapply(gr_list,"[[","code"))
}


#plot how best performance varies as a line gets some links substituted by rings
# diag is the default optimization method since it ise the only time-feasable one for many dimensions
plot_line_to_CHchain = function(L, ring_size= 3,space = 1, n_opt_mode = "diag")
{

  #number of link the subsitution unit will occupy
  n_sub_link = ring_size%/%2

  #create substitution unit
  unit_str = "Lx+Cy"
  unit_str = gsub("x",space,gsub("y",ring_size,unit_str))
  unit_gr = gfc(unit_str)

  print(unit_str)

  gr_list = list()

  gr_list[[1]] = line(L)

  #subsitute 1 to L link to unit
  for(i in 1:L)
  {
    size = i*n_sub_link + (i-1)*space
    if(size> (L-1))
      break

    sx = (L+1-size)%/%2 + (L+1-size)%%2
    dx = (L+1-size)%/%2

    #print(paste(sx,dx,size))

    #create substitution unit
    subst = chain( ring(ring_size), rep = i ,space,HANDLES = F)

        #add left&right appendix
    temp_gr = join_nolink(line(sx),subst)
    temp_gr = join_nolink(temp_gr,line(dx))

    gr_list[[i+1]] = temp_gr

  }

  N = length(gr_list)

  plot(-10,-10,ylim = c(0,1.4), xlim= c(0,N),xlab = "subbed unit", ylab=expression(p[max]))
	legend("topleft", title = paste("L",L," w\ C", ring_size, " ring subst",sep=""),
         legend = c("gr", " '' nophi"),
         col = c(2,2),lty= c(1,3), pch = c(19,-1))
  grid()

  best = opt_and_plot_list(gr_list,opt_mode = n_opt_mode,mode = "range", range = 0:(N-1),nophi = TRUE, col = 2, pch = 19)

  print( paste("Massimo di trasporto:",best))

  return(sapply(gr_list,"[[","code"))

}

#see how a handled chain of unit (glued by add_nolink) scales
plot_unit_scaling = function(gr_unit, range = 1:10, space = 0, HANDLES = TRUE, n_opt_mode = "diag", n_mode = "sites", ADD = FALSE)
{

  gr_list = lapply(range, chain,gr = gr_unit, space = space, HANDLES = HANDLES)

  lims = get_list_range(gr_list,n_mode)

  if(!ADD)
  {
    plot(-10,-10,ylim = c(0,1.4), xlim = lims,xlab = n_mode, ylab="p_max")
    legend("topright", title = paste(gr_unit$code," chain scaling",sep=""),
           legend = c(gr_unit$code, " '' nophi"),
           col = c(2,2),lty= c(1,3), pch = c(19,-1))
    grid()
  }


  best = opt_and_plot_list(gr_list,opt_mode = n_opt_mode ,mode = n_mode, range = range,nophi = TRUE, col = 2, pch = 19)

  print( paste("Massimo di trasporto:", best))

}

#see how a handled chain of unit (glued by add_nolink) scales with time
plot_unit_t_scaling = function(gr_unit, range = 1:10, space = 0, HANDLES = TRUE, n_opt_mode = "diag", n_mode = "sites", L_LINE = TRUE)
{

  gr_list = lapply(range, chain,gr = gr_unit, space = space, HANDLES = HANDLES)

  lims = get_list_range(gr_list,n_mode)
  t_max = max_search_time(tail(gr_list,1)[[1]])

  t_progression_plot(ymax = t_max,xmax = lims[2], x_lab = n_mode)
	legend("topleft", title = paste(gr_unit$code," chain t-scaling",sep=""),
         legend = c(gr_unit$code, " '' nophi", "L"),
         col = c(2,2,3),lty= c(1,3,1), pch = c(19,-1,-1))

  best = opt_and_t_plot_list(gr_list,opt_mode = n_opt_mode ,mode = n_mode, range = range,nophi = TRUE, col = 2, pch = 19)

  if(L_LINE)
    abline(line_t_lm(),col = 3)



  print( best)

}

#plot how just the size of the cycle unit influence the transport performance
plot_unit_size_scaling = function( range = 3:10, rep = 10, space =0, HANDLES = TRUE,  n_opt_mode = "diag", n_mode = "dist", L_LINE = TRUE)
{
  gr_base = lapply(range,ring)
  gr_list = lapply(gr_base, chain,rep = rep, space = space, HANDLES = HANDLES)

  lims = get_list_range(gr_list,n_mode)

  progression_plot(paste("cycle size", rep, " unit chain scaling",sep=""),xmax = lims[2], x_lab = n_mode)
	legend("topleft",
         legend = c("chain size prog", "L"),
         col = c(2,3),lty= c(-1,1), pch = c(19,-1))

  best = opt_and_plot_list(gr_list,opt_mode = n_opt_mode ,mode = n_mode, range = range,nophi = F, LINES = F, col = 2, pch = 19)

  if(L_LINE)
    plot_line_trend(lims[1]:lims[2])



  print( best)

}

plot_chain_freature_comp = function(unit, range = 1:10, L_LINE = FALSE)
{
  gr_list_best =  lapply(range, chain,gr = unit, space = 0, HANDLES = TRUE)
  gr_list_space = lapply(range, chain,gr = unit, space = 1, HANDLES = TRUE)
  gr_list_noH =   lapply(range, chain,gr = unit, space = 0, HANDLES = FALSE)
  gr_list_worst = lapply(range, chain,gr = unit, space = 1, HANDLES = FALSE)

  lims = get_list_range(gr_list_worst,"dist")

  progression_plot(xmax = lims[2], x_lab = "dist")
	legend("topleft",
         legend = c("best","1 space", "no H", "1 sp, noH"), ncol = 2,
         col = 2:5,lty= 1, pch = 19)

  best = opt_and_plot_list(gr_list_best, opt_mode = "diag" ,mode = "dist",nophi = F, LINES = T, col = 2, pch = 19)
  opt_and_plot_list(gr_list_space,       opt_mode = "diag" ,mode = "dist",nophi = F, LINES = T, col = 3, pch = 19)
  opt_and_plot_list(gr_list_noH,         opt_mode = "diag" ,mode = "dist",nophi = F, LINES = T, col = 4, pch = 19)
  opt_and_plot_list(gr_list_worst,       opt_mode = "diag" ,mode = "dist",nophi = F, LINES = T, col = 5, pch = 19)

  if(L_LINE)
    plot_line_trend(lims[1]:lims[2])



  print( best)
}

plot_chain_size_scaling = function(range = 3:10 ,size = 10, HANDLES = TRUE, space = 0, L_LINE = TRUE)
{
  gr_base = lapply(range,ring)

  gr_list = lapply(gr_base, chain, rep = size, HANDLES = HANDLES, space = space)
  lims = get_list_range(gr_list,"dist")

  progression_plot(xmax = max(range), x_lab = "unit size")
	legend("topleft",
         legend = c( paste("C chain",size,"u"), "L"),
         col = 2:3,lty= c(-1,1), pch = c(19,-1) )


  opt_and_plot_list(gr_list, opt_mode = "diag" ,mode = "range", range = range, nophi = F, LINES = F, col = 2, pch = 19)

  if(L_LINE )
  {
    #plot same lenght line graph performance
    dist_vec = sapply(gr_list, distance)
    line_list = lapply(dist_vec, line)

    p_vec = sapply( line_list, p_max_func)

    lines(range,p_vec, col = 3)

  }
}

plot_Rc_size_scaling = function(range = (2:6) *2,size = 10, HANDLES = TRUE, space = 0, L_LINE = TRUE)
{

  gr_base = lapply(range, correct_cut_square_unit)

  gr_list = lapply(gr_base, chain, rep = size, HANDLES = HANDLES, space = space)
  lims = get_list_range(gr_list,"dist")

  progression_plot(xmax = max(range)+2, x_lab = "N")
	legend("topleft",
         legend = c(paste("RcN chain",size,"u"), "L"),
         col = 2:3,lty= c(-1,1), pch = c(19,-1)  )

  opt_and_plot_list(gr_list, opt_mode = "diag" ,mode = "range", range = range, nophi = F, LINES = F, col = 2, pch = 19)

  if(L_LINE)
    abline( h = p_max(line(size+3)), col = 3)

}

plot_ch_vs_couples = function( d_max = 30, L_LINE = TRUE )
{

  gr_list_ch = lapply(3:((d_max-2)*2), ringNhandle)

  gr_base = lapply(3:(d_max-2), ring)
  gr_list_cch = lapply(gr_base, chain, rep = 2, HANDLES = TRUE,space = 0)

  lims = get_list_range(gr_list_cch,"dist")

  progression_plot(xmax = d_max, x_lab = "dist")
	legend("topright",
         legend = c("Ch", "L1+Ca/Ca+L1", "L"),
         col = c(2,4,3),lty= 1, pch = c(19,19,-1)  )

  best =  opt_and_plot_list(gr_list_ch , mode = "dist",nophi = F, LINES = T, col = 2, pch = 19)
          opt_and_plot_list(gr_list_cch, mode = "dist",nophi = F, LINES = T, col = 4, pch = 19)

  if(L_LINE)
    plot_line_trend(3:d_max)
}

plot_chain_vs = function( d_max = 30, L_LINE = TRUE )
{

  gr_list_c3 =  lapply(1:(d_max-2),     chain, gr = ring(3),                   HANDLES = TRUE, space = 0)
  gr_list_cc4 = lapply(1:(d_max-2),     chain, gr = correct_cut_square_unit(), HANDLES = TRUE, space = 0)
  gr_list_c4 =  lapply(1:((d_max-2)/2), chain, gr = ring(4),                   HANDLES = TRUE, space = 0)

  progression_plot(xmax = d_max, x_lab = "dist")
	legend("topright",
         legend = c("C3 u", "C4 u", "Rc4(1-3) u", "L"), ncol = 2,
         col = c(2,4,5,3),lty= 1 )

  opt_and_plot_list(gr_list_c3,  opt_mode = "diag", mode = "dist",nophi = F, POINTS = F, LINES = T, col = 2)
  opt_and_plot_list(gr_list_c4,  opt_mode = "diag", mode = "dist",nophi = F, POINTS = F, LINES = T, col = 4)
  opt_and_plot_list(gr_list_cc4, opt_mode = "diag", mode = "dist",nophi = F, POINTS = F, LINES = T, col = 5)

  if(L_LINE)
    plot_line_trend(3:d_max)
}

plot_chain_vs_ch = function( d_max = 30, L_LINE = TRUE )
{

  gr_list_c3 =  lapply(1:(d_max-2),     chain, gr = ring(3),                   HANDLES = TRUE, space = 0)
  gr_list_ch =  lapply(3:((d_max-2)*2),   ringNhandle)

  progression_plot(xmax = d_max, x_lab = "dist")
	legend("topright",
         legend = c("C3 u", "Ch","L"), ncol = 1,
         col = c(2,4,3),lty= 1, pch = c(19,19,-1) )

  opt_and_plot_list(gr_list_c3,  opt_mode = "diag", mode = "dist",nophi = F, POINTS = T, LINES = T, col = 2, pch = 19)
  opt_and_plot_list(gr_list_ch,  opt_mode = "",     mode = "dist",nophi = F, POINTS = T, LINES = T, col = 4, pch = 19)

  if(L_LINE)
    plot_line_trend(3:d_max)
}

