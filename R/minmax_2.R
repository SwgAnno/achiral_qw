# 2 dimensional optimization methods

# still depends upon 1D opti in minmax.R


#abstraction impose a special method for 2 free parameter performance
performance_2 = function(no_phi_gr, sample = seq(0,2*pi,.5))
{
  if( length( no_phi_gr$re_coord)!=2 )
  {
    print("Numero di fasi libere non corretto")
    return(-1)
  }

  n = length(sample)
  phi_vec = exp(sample *1i)

  phi_mat = expand.grid(phi_vec,phi_vec)

  gr_vec = apply(phi_mat,1,rephase,gr=no_phi_gr)

  ret = sapply(gr_vec,p_max_func)

	dim(ret) = c(n,n)

	return(ret)
}

performance_diag = function(no_phi_gr, sample = seq(0,2*pi,.1) )
{

  n = length(no_phi_gr$re_coord)
  phi_vec = exp(sample *1i)

  phi_mat = matrix( rep(phi_vec,n), ncol = n)

  gr_vec = apply(phi_mat,1,rephase,gr=no_phi_gr)

  ret = sapply(gr_vec,p_max_func)

	return(ret)
}

t_performance_diag = function(no_phi_gr, sample = seq(0,2*pi,.1) )
{

  n = length(no_phi_gr$re_coord)
  phi_vec = exp(sample *1i)

  phi_mat = matrix( rep(phi_vec,n), ncol = n)

  gr_vec = apply(phi_mat,1,rephase,gr=no_phi_gr)

  ret = sapply(gr_vec, t_p_max_func)

	return(ret)
}

t_performance_2 = function(no_phi_gr, sample = seq(0,2*pi,.1))
{
  if( length( no_phi_gr$re_coord)!=2 )
  {
    print("Numero di fasi libere non corretto")
    return(-1)
  }

  n = length(sample)
  phi_vec = exp(sample *1i)

  phi_mat = list()

  phi_mat = expand.grid(phi_vec,phi_vec)

  gr_vec = apply(phi_mat,1,rephase,gr=no_phi_gr)

  ret = sapply(gr_vec,t_p_max_func)

	dim(ret) = c(n,n)

	return(ret)
}




plot_performance_2 = function(gr, by = .1)
{
  prog = seq(0,2*pi,by)

  p_mat = performance_2(gr, prog)

  if(length(p_mat) == -1)
  {
    print("Errore fasi libre plot_performance_2")
    return (-1)
  }

  colnames(p_mat) = prog
  rownames(p_mat) = prog

  print( paste( "Massimo di trasporto", max(p_mat)))

  pos_max = which.max(p_mat) -1
  l = length(prog)
  phi1 = prog[pos_max%%l +1]
  phi2 = prog[pos_max%/%l +1]
  print( paste( "con fasi", phi1, phi2))

  #this way od doing things keeps phi1 on the x

  image.plot(prog, prog, p_mat, xlab = expression(phi[1]), ylab = expression(phi[2]))

}

plot_performance_diag = function(gr, by = .1, t_col = FALSE)
{

  prog = seq(0, 2*pi, by)

  max_vec = performance_diag( gr, prog)

  print( paste( "Massimo di trasporto", max(max_vec)))

  pos_max = which.max(max_vec)
  print( paste( "Con fase: ", prog[pos_max]))

  if(t_col)
  {
    t_vec   = t_performance_diag( gr, prog)

    #values for plot legend
    col_scale_leg = round( seq(min(t_vec), max(t_vec),length.out = 8) , digits = 1)

    #project t_vec onto 0-1 from min to max
    t_vec = (t_vec - min(t_vec))/ (max(t_vec) - min(t_vec))

    col_ramp = colorRamp(c( "red", "blue"))
    col_vec = rgb( col_ramp(t_vec)/255 )

    performance_plot()
    legend( "topright",
            legend = col_scale_leg,
            col = rgb( col_ramp(seq(0,1,.1))/255 ),
            pch = 20,
            text.width = .5,
            ncol = 2)

    points(prog, max_vec, col = col_vec, pch = 20)
  }

  else
  {
    performance_plot()

    legend( "topright",
            legend = gr$code,
            col = 2,
            lty = 1,
            )
    lines(prog, max_vec, col = 2,lty = 1)
  }

}

plot_t_performance_2 = function(gr, by = .1, SURFACE = FALSE)
{
  prog = seq(0,2*pi,by)

  p_mat = t_performance_2(gr, prog)

  if(length(p_mat) == -1)
  {
    print("Errore fasi libre plot_performance_2")
    return (-1)
  }

  colnames(p_mat) = prog
  rownames(p_mat) = prog

  if(SURFACE)
    rotate.persp(prog, prog, p_mat)
  else
    image.plot(prog, prog, p_mat, xlab = expression(phi[1]), ylab = expression(phi[2]))

}

debug_opti = function( gr , by = .1, ADD_GRAPH = FALSE)
{
  n = length( gr$re_coord)

  if( !ADD_GRAPH)
    plot_performance(gr,by)

  get_perf = function(phi)
    p_max_func( rephase(gr, exp(1i*phi)  ))

  init = matrix( runif(GLOBAL_RANDOM_START_COUNT*n) ,ncol = n)
  init = init * 2*pi

  #optim control parameters
  my_control = list()
    my_control$fnscale = -1


  max = apply(init,1 ,optim, get_perf, method = "L-BFGS-B", lower = 0,upper = 2*pi, control = my_control)

  max_phi = sapply( max, "[[", "par")
  print(max_phi)

  max_p = apply(max_phi,2,get_perf)
  print(max_p)

  points( max_phi[1,],max_phi[2,], col = "white", pch = 20)
  points( init[,1], init[,2], col = "black", pch= 20)

  print( paste( "Massimo algoritmo: ", max(max_p)))
  print( paste( "Con fasi:",  max_phi[,which.max(max_p)])  )

}

plot_cut_grid =function(N = 8, HANDLES = FALSE, SURFACE = FALSE)
{

  prog= 1:N
  cut_vec = expand.grid(prog,prog)

  gr_list = apply(cut_vec, 1, ringNcut, N = N, HANDLES = HANDLES)
  gr_list = lapply(gr_list, optimize_gr)

  p_mat = sapply(gr_list, p_max_func)

  dim(p_mat) = c(N,N)
  colnames(p_mat) = prog
  rownames(p_mat) = prog

  if(SURFACE)
    rotate.persp(prog, prog, p_mat)
  else
    image.plot(prog, prog, p_mat, xlab = "cut start", ylab = "cut end")
}

plot_t_cut_grid =function(N = 8, HANDLES = FALSE, SURFACE = FALSE)
{

  prog= 1:N
  cut_vec = expand.grid(prog,prog)
  gr_list = lapply(gr_list, optimize_gr)

  gr_list = apply(cut_vec, 1, ringNcut, N = N, HANDLES = HANDLES)

  t_mat = sapply(gr_list, t_p_max_func)

  dim(t_mat) = c(N,N)
  colnames(t_mat) = prog
  rownames(t_mat) = prog

  if(SURFACE)
    rotate.persp(prog, prog, t_mat)
  else
    image.plot(prog, prog, t_mat, xlab = "cut start", ylab = "cut end")
}

plot_topo_grid = function(op = "/",HANDLES = TRUE, range = 1:10,target = "p" ,SURFACE= TRUE, opt_mode ="")
{
  base_str = paste("Ca" ,op,"Cb", sep = "")

  if(HANDLES)
  {
    base_str = paste("L1+", base_str, "+L1", sep = "")
  }

  max_func = function(gr)
       	return(p_max_func(optimize_gr(gr, mode = opt_mode)))

  if(target == "t")
  {
    max_func = function(gr)
          return(t_p_max_func(optimize_gr(gr,mode = opt_mode)))
  }


  # support function to extract pairwise best performance
  pair_perf = function( pos)
  {
    gr = gfc (gsub("a",pos[1],gsub("b",pos[2],base_str)))

    return(max_func(gr))
  }

  l = length(range)
  mat = expand.grid(range,range)

  p_mat = apply(mat,1,pair_perf)
  dim(p_mat) = c( l, l)

  colnames(p_mat) = range
  rownames(p_mat) = range

  print(p_mat)

  if(SURFACE)
  {
        col.pal<-colorRampPalette(c("blue", "red"))
    colors<-col.pal(100)
    # height of facets
    z.facet.center <- (p_mat[-1, -1] + p_mat[-1, -ncol(p_mat)] + p_mat[-nrow(p_mat), -1] + p_mat[-nrow(p_mat), -ncol(p_mat)])/4
    # Range of the facet center on a 100-scale (number of colors)
    z.facet.range<-cut(z.facet.center, 100)

    persp(range, range, p_mat,
      theta = -60, phi = 30,
      xlab = "ring 1",ylab = "ring 2", zlab = "t p_max",
      shade=NA, col=colors[z.facet.range], border="grey80",
      box=TRUE, axes = TRUE, nticks = 5, ticktype = "detailed")
  }
    #rotate.persp(range, range, p_mat)
  else
    image.plot(range, range, p_mat, xlab = "phi_1", ylab = "phi_2")

}

#plot the scaling of a couple with costant sum of ring sites
plot_isobare_progression = function( N, glue = "/", HANDLES = TRUE)
{
  base_str = paste("Ca" ,glue,"Cb", sep = "")

  if(HANDLES)
  {
    base_str = paste("L1+", base_str, "+L1", sep = "")
  }

  #print(base_str)

  n_mat = matrix( c(1:(N-1), (N-1):1), ncol = 2)

  ab_sub = function(subs, str)
    gsub("a",subs[1],gsub("b",subs[2],str))

  str_vec = apply(n_mat,1,ab_sub, str = base_str)

  print(str_vec)

  gr_list = lapply( str_vec, gfc)

  plot(-10,-10,ylim = c(0,1.4), xlim = c(1,(N-1)),xlab = "first ring sites", ylab="p_max")
	legend("topright", title = paste(N, "sum"),
         legend = c("isobare couple perf.", " '' nophi"),
         col = c(2,2),lty= c(1,3), pch = c(19,-1))
  grid()

  opt_and_plot_list(gr_list,mode = "range", range = 1:(N-1),nophi = TRUE, col = 2, pch = 19)


}

plot_N_couple_progression = function( range = 2:10, n_mode = "dist")
{

  if(n_mode=="dist")
    label = "min_dist"
  else
    label = "N"

  gr_base = lapply(range,ring)

  ######
  # add couples

  gr_list = lapply(gr_base,chain, rep = 2, space = 1, HANDLES = TRUE)
  xmax = max(get_list_range(gr_list, n_mode))
  {
  plot(-10,-10,ylim = c(0,1.3), xlim= c(1, xmax),xlab = n_mode, ylab=expression(p[max]))
	legend("topright", legend = c("CN+CN"," '' no handles", "CN/CN"," '' no handles" ), col = c(2,2,3,3), lty = 1:2, pch = 1)
  grid()
  }

  opt_and_plot_list(gr_list,mode = n_mode, col = 2)

  gr_list = lapply(gr_base,chain, rep = 2, space = 1, HANDLES = FALSE)
  opt_and_plot_list(gr_list,mode = n_mode, col = 2, lty = 2)

  ######
  # glue couples

  gr_list = lapply(gr_base,chain, rep = 2, space = 0, HANDLES = TRUE)
  opt_and_plot_list(gr_list,mode = n_mode, col = 3)

  gr_list = lapply(gr_base,chain, rep = 2, space = 0, HANDLES = FALSE)
  opt_and_plot_list(gr_list,mode = n_mode, col = 3, lty = 2)

}

plot_t_N_couple_progression = function( range = 2:10, n_mode = "dist")
{

  if(n_mode=="dist")
    label = "min_dist"
  else
    label = "N"

  gr_base = lapply(range,ring)

  ######
  # add couples

  gr_list = lapply(gr_base,chain, rep = 2, space = 1, HANDLES = TRUE)
  xmax = max(get_list_range(gr_list, n_mode))
  tmax = max_search_time(tail(gr_list,1)[[1]])
  {
  plot(-10,-10,ylim = c(0,tmax), xlim= c(1, xmax),xlab = label, ylab="t")
	legend("topleft", legend = c("CN+CN"," '' /w handles", "CN/CN"," '' /w handles" ), col = c(2,2,3,3), lty = 1:2)
  grid()
  }

  opt_and_t_plot_list(gr_list,mode = n_mode, col = 2)

  gr_list = lapply(gr_base,chain, rep = 2, space = 1, HANDLES = FALSE)
  opt_and_t_plot_list(gr_list,mode = n_mode, col = 2, lty = 2)

  ######
  # glue couples

  gr_list = lapply(gr_base,chain, rep = 2, space = 0, HANDLES = TRUE)
  opt_and_t_plot_list(gr_list,mode = n_mode, col = 3)

  gr_list = lapply(gr_base,chain, rep = 2, space = 0, HANDLES = FALSE)
  opt_and_t_plot_list(gr_list,mode = n_mode, col = 3, lty = 2)

}

plot_couple_space_progression = function(range = 0:10,base = ring(4), HANDLES = FALSE)
{
  gr_list = lapply(range,chain,gr = base,rep = 2,space = 0,HANDLES = HANDLES)
  list_range = get_list_range(gr_list, "dist")



  progression_plot("", xmax = list_range[2],x_lab = "dist")

  opt_and_plot_list(gr_list, mode = "dist",col = 2)
  plot_line_trend(list_range[1]:list_range[2] +1)
}

plot_ringNc_short_prog = function(range = 3:10, HANDLES = FALSE)
{
  t_progression_plot("",ymax = max(range), xmax = max(range),x_lab = "ring size")
  legend("topleft",legend = c("Rc FM","Rc TC1"), col = 2:3, pch = 19, lty = 1)

  pushFC(T)

  gr_list = lapply(range, ringNcut , HANDLES = F)
  opt_and_t_plot_list(gr_list, range = range,col = 2, pch = 19)


  popFC()
  pushFC(F)
  pushTC(1)

  opt_and_t_plot_list(gr_list, range = range,col = 3, pch = 19)

  popFC()
  popTC()
}
