#gets back time stamp of maxima and minima
deriv_roots = function(gr, from = 0,to = max_search_time(gr))
{

	f = function(t) deriv_p_l_t(gr,t)

	out = uniroot.all(f,lower = from, upper = to, tol = GLOBAL_UNIROOT_TOL)
	if(length(out)== 0)
  {
    #print("WARNING: no maxima found")
    return(to)
  }
	return( out  )


}

# stone age stationary point search (aka poitwise search of derivative change of sign)
# minmax arg can be used to get only max or mins if it's accordingly >0 or <0
deriv_roots_uga_buga = function(gr, from = 0,to = max_search_time(gr), by = .001, minmax = 0)
{
  t_seq = seq(from,to,by)

	evo = deriv_full(gr,from,to,by)

	#shift evo right and add dummy
	p1 = c(2,head(evo,-1))

  if(minmax > 0)
    stationary =  evo>=0 & p1<0

  else if(minmax < 0)
    stationary =  evo<=0 & p1>0

  else
    stationary = evo*p1 < 0 | evo == 0

  return( t_seq[stationary])

}

p_max = function(gr, from = 0,to = max_search_time(gr))
{
  if(YOLO)
    candidates = deriv_roots_uga_buga(gr,from,to)
  else
    candidates = deriv_roots(gr,from,to)

	p_vec = Mod(evo_l_t(gr,gr$start, candidates))^2

	return( max(p_vec[gr$target,]))
}

p_max_first = function(gr, from = distance(gr)/4, to = distance(gr)*2)
{
  if(YOLO)
    candidates = deriv_roots_uga_buga(gr,from,to)
  else
    candidates = deriv_roots(gr,from,to)

	p_vec = Mod(evo_l_t(gr,gr$start, candidates))^2

  return(p_vec[gr$target,1])
}

p_lm = function(gr ,model = get_topo_lm(gr))
{
  p_vec = Mod( evo_l_lm(gr, gr$start,model) )^2

  return(p_vec[gr$target])
}

t_p_max = function(gr, from = 0,to = max_search_time(gr))
{
  if(YOLO)
    candidates = deriv_roots_uga_buga(gr,from,to,minmax = 1)
  else
    candidates = deriv_roots(gr,from,to)

	p_vec = Mod(evo_l_t(gr,gr$start, candidates))^2

	pos_max = which.max(p_vec[gr$target,])
	return( candidates[pos_max])
}

t_p_max_first = function(gr,  from = distance(gr)/4, to = distance(gr)*3/4)
{
  if(YOLO)
    candidates = deriv_roots_uga_buga(gr,from,to,minmax = 1)
  else
    candidates = deriv_roots(gr,from,to)

	return( candidates[1])
}


# do both p_max and t_p_max
locate_max = function(gr, from = 0,to = max_search_time(gr))
{
  candidates = deriv_roots(gr,from,to)
	p_vec = Mod(evo_l_t(gr,gr$start, candidates))^2

	pos_max = which.max(p_vec[gr$target,])

  out = list()
	out$m = max(p_vec[gr$target,])
	out$t = candidates[pos_max]
	return(out )
}

evo_vs_deriv = function(gr,from = 0, to= max_search_time(gr), by = .1)
{
	t_seq = seq(from,to,by)

	p_vec = Mod(evo_full(gr, from, to ,by))^2
	deriv_vec = deriv_p_l_t(gr,t = t_seq)

	plot(0,0,ylim = c(-2,2), xlim= c(from,to ),xlab = "l", ylab="p(l)")
  grid()

	lines(t_seq, p_vec[gr$target,]/sum(p_vec[,1]), col= gr$target)
	lines(t_seq, deriv_vec,col = gr$target+1)

	legend(3,1, legend = c("P(t)","P'(t)"), col = c(gr$target,gr$target+1)  ,lty = 1)
}

performance = function( no_phi_gr,sample = seq(0,2*pi,.1))
{
	phi_vec = exp(sample *1i)

	gr_vec = lapply(phi_vec,rephase,gr=no_phi_gr)

	ret = sapply(gr_vec, p_max_func)

	return(ret)
}

lm_performance = function(no_phi_gr,sample = seq(0,2*pi,.1), lm = get_topo_lm(no_phi_gr))
{
  phi_vec = exp(sample *1i)

	gr_vec = lapply(phi_vec,rephase,gr=no_phi_gr)

	ret = sapply(gr_vec,p_lm,lm)

  return(ret)
}

t_performance = function(no_phi_gr, sample = seq(0,2*pi,.1))
{
	phi_vec = exp(sample *1i)

	gr_vec = lapply(phi_vec,rephase,gr=no_phi_gr)

  ret = sapply(gr_vec,t_p_max_func)

	return(ret)
}

plot_performance = function(gr, by = .1, t_col = FALSE)
{
  l = length(gr$re_coord)
  if(l==1)
    plot_performance_1( gr, by, t_col)
  else if(l==2)
    plot_performance_2( gr, by)
  else
    print("plot_performance not supported for target phase DOF")
}

#plot probability performance along with "time performance" information displayed with color

plot_performance_1 = function(gr, by = .1, t_col = FALSE)
{
  prog = seq(0, 2*pi, by)

  max_vec = performance( gr, prog)

  if(t_col)
  {
    t_vec   = t_performance( gr, prog)

    #values for plot legend
    col_scale_leg = format( seq(min(t_vec), max(t_vec),length.out = 10) , digits = 2)

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

#plots the performance of all the graphs in the gr list
plot_performance_multi_1 = function( gr_list, by = .1)
{
  #convert to one element list if you get just a graph
  if(class(gr_list)== "graph")
    gr_list = list(gr_list)

	prog = seq(0,2*pi,by)

	performance_plot()

	codes = sapply(gr_list, "[[", "code")
	legend(.2,1.3, legend = codes, col = 1:length(gr_list), pch = 1)

	max_vec = lapply(gr_list,performance ,sample = prog)

	for(i in seq_along(max_vec))
    lines( prog, max_vec[[i]], col = i)

}

plot_t_performance = function(gr, by = .1)
{
  l = length(gr$re_coord)
  if(l==1)
    plot_t_performance_1( gr, by )
  else if(l==2)
    plot_t_performance_2( gr, by)
  else
    print("plot_t_performance not supported for target phase DOF")
}

#phase depenent time-of-arrival plot
plot_t_performance_1 = function( gr, by = .1)
{
  #convert to one element list if you get just a graph
  if(class(gr)== "graph")
    gr = list(gr)

	prog = seq(0,2*pi,by)
	t_max = max(sapply(gr,max_search_time))

	t_performance_plot(t_max = t_max)
	codes = sapply(gr, "[[", "code")
	legend(.2,t_max, legend = codes, col = 1:length(gr), pch = 1)

	max_vec = lapply(gr,t_performance ,prog)

	for(i in seq_along(max_vec))
    lines( prog, max_vec[[i]], col = i)

}

#DEPRECATED: now the sam thing is done better with first maximum search

plot_performance_vs_lm = function(gr, by = .1)
{

  prog = seq(0,2*pi,by)

	performance_plot()

	lg = c( paste(gr$code, " lm"),paste(gr$code, " best"))
	legend("topleft", legend = lg, col = 1:2, pch = 1)

	####
	max_vec = lm_performance(gr,prog)
	lines( prog, max_vec, col = 1)

  ####
  max_vec = performance(gr,prog)
  lines( prog, max_vec, col = 2)

}

#overall best maximum performance vs firt maximum performance

plot_performance_TC_vs_fist_max = function(gr,TC_vec = getTC(), by = .1)
{

  prog = seq(0,2*pi,by)

	performance_plot()

	lg = c( paste(gr$code, " first"),paste(gr$code, " best, TC", TC_vec))
	legend("topleft", legend = lg, col = 1:(1+ length(TC_vec)), pch = 1)

	#"push" new value
  pushFC(TRUE)

	####
	max_vec = performance(gr,prog)
	lines( prog, max_vec, col = 1)

  setFC(FALSE)

  for(i in seq_along(TC_vec))
  {
    pushTC(TC_vec[i])

    max_vec = performance(gr,prog)
    lines( prog, max_vec, col = i+1)

    popTC()

  }

  popFC()

}

plot_even_progression = function(num = 5)
{
	prog = seq(0,2*pi,.1)
	nums = seq(2,2*num,2)

	plot(-10,-10,ylim = c(0,1.3), xlim= c(0,2*pi),xlab = "e^ phi / (0 -> 2pi)", ylab="p_max")
  grid()

	labels = sapply(nums,nsub,"C(N)")
	legend(.2,1.3, legend = labels, col = nums, pch = 1)

	for( i in seq_along(nums))
	{
		gr = ring(nums[i])
		max_vec = performance(gr, sample = prog)

		lines( prog, max_vec, col = nums[i]/2)
	}

}

plot_vs_ring = function(N, by = .1)
{
	prog = seq(0,2*pi,by)

	plot(-10,-10,ylim = c(0,1.3), xlim= c(0,2*pi),xlab = "e^ phi / (0 -> 2pi)", ylab="p_max")
  grid()

	labels = nsub(N, c("C(N)","Ch(N)"))
	legend(.2,1.3, legend = labels, col = 2:3, pch = 1)

	gr = ring(N)
	max_vec = performance(gr, sample = prog)
	lines( prog, max_vec, col = 2)

	gr = ringNhandle(N)
	max_vec = performance(gr, sample = prog)
	lines( prog, max_vec, col = 3)

}

plot_N_progression = function(g_max  = 10, mode = "sites")
{
	prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,1.3), xlim= c(0,g_max),xlab = mode, ylab = expression(p[max]))
	legend("topright", legend = c("C", "Ch","L"), col = 1:3, pch = 19)
  grid()

  ####
  if(mode == "dist")
    v_max = g_max*2
  else
    v_max = g_max

	gr_vec = lapply(3:v_max,ring)
	opt_and_plot_list(gr_vec,mode =mode ,nophi = F, col = 1, pch = 19)

  ####
	if(mode == "dist")
    v_max = (g_max-2)*2
  else
    v_max = g_max-2

	gr_vec = lapply(3:v_max,ringNhandle)
	opt_and_plot_list(gr_vec,mode =mode ,nophi = F, col = 2, pch = 19)


	####
  if(mode == "dist")
    v_max = g_max+1
  else
    v_max = g_max

	gr_vec = lapply(2:v_max,line)
	opt_and_plot_list(gr_vec, opt_mode = "none",mode =mode ,nophi = F, col = 3, pch = 19)



}

plot_t_N_progression = function(g_max  = 10, mode = "sites")
{
	prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,g_max*getTC()), xlim= c(0,g_max),xlab = mode, ylab= expression(t[p[max]]))
	legend("topleft", legend = c("C", "Ch","L"), col = 1:3, pch = 19)
  grid()


	####
	if(mode == "dist")
    v_max = g_max*2
  else
    v_max = g_max

	gr_vec = lapply(3:v_max,ring)

	opt_and_t_plot_list(gr_vec,mode =mode ,nophi = F, col = 1, pch = 19)

	####
	if(mode == "dist")
    v_max = (g_max-2)*2
  else
    v_max = g_max-2

	gr_vec = lapply(3:v_max,ringNhandle)

	opt_and_t_plot_list(gr_vec,mode =mode ,nophi = F, col = 2, pch = 19)

	####
  if(mode == "dist")
    v_max = g_max+1
  else
    v_max = g_max

	gr_vec = lapply(2:v_max,line)

	opt_and_t_plot_list(gr_vec,opt_mode = "none", mode =mode ,nophi = F, col = 3, pch = 19)

}

plot_topo_progression = function(code = "C",range = 2:15)
{
  prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,1.3), xlim= range(range),xlab = "N", ylab="p_max_phi")
	legend("topleft", legend = paste(code,"N",sep = ""), col = 2, pch = 1)
  grid()

  gr_vec = get_N_topo_opt(code ,range)

  max_vec = sapply(gr_vec, p_max_func)

	lines( range, max_vec, pch = 19, col = 1)
	points( range, max_vec, pch = 19, col = 1)


}

# increasing size graph of the selected topo
# comparison between phase-optimized and not optimize graph
plot_topo_vs_nophase = function(code = "C", range = 2:15)
{
  prog = seq(0,2*pi,.1)

  progression_plot(xmax = max(range))

  lg = c( code,paste(code," opt'd",sep = ""))
	legend("topleft", legend = lg, col = 2:3, pch = 19, lty = 1)

  gr_vec = get_N_topo(code ,range)

  ####
  max_vec = sapply(gr_vec, p_max_func)

	lines( range, max_vec, pch = 19, col = 2)
	points( range, max_vec, pch = 19, col = 2)

	####
  gr_vec = lapply(gr_vec,optimize_gr)

  max_vec = sapply(gr_vec, p_max_func)

  lines( range, max_vec, pch = 19, col = 3)
	points( range, max_vec, pch = 19, col = 3)

}

# increasing size graph of the selected topo
#comparison between overall best maxima and first maximum performance
plot_topo_vs_first_max = function(code = "C", range = 2:15)
{
  prog = seq(0,2*pi,.1)


  lg = c( paste(code,"N first",sep = ""),paste(code,"N best",sep = ""))
	legend("topleft", legend = lg, col = 2:3, pch = 1)

  performance_plot()

  gr_vec = get_N_topo(code ,range)



  ####
  pushFC(TRUE)
  gr_vec = lapply(gr_vec,optimize_gr)

	max_vec = sapply(gr_vec, p_max_func)

	lines( range, max_vec, pch = 19, col = 2)
	points( range, max_vec, pch = 19, col = 2)

	####
	setFC(FALSE)
  gr_vec = lapply(gr_vec,optimize_gr)

	max_vec = sapply(gr_vec, p_max_func)

  lines( range, max_vec, pch = 19, col = 3)
	points( range, max_vec, pch = 19, col = 3)


  popFC()

}

#DEPRECATED: now the sam thing is done better with first maximum search

# increasing size graph of the selected topo
#comparison between linear model maxima and overall best maxima
plot_topo_vs_lm = function(code = "C", range = 2:15)
{
  prog = seq(0,2*pi,.1)


  lg = c( paste(code,"N",sep = ""),paste(code,"N lm",sep = ""))
	legend("topleft", legend = lg, col = 2:3, pch = 1)

  performance_plot()

  gr_vec = get_N_topo(code ,range)

  ####
  gr_vec = lapply(gr_vec,optimize_gr)

  max_vec = sapply(gr_vec, p_max_func)

	lines( range, max_vec, pch = 19, col = 2)
	points( range, max_vec, pch = 19, col = 2)

	####

	model = get_topo_lm(gr_vec[[1]])
  gr_vec = lapply(gr_vec,optimize_gr, model = model, mode = "lm")
  max_vec = sapply(gr_vec, p_lm,model)

  lines( range, max_vec, pch = 19, col = 3)
	points( range, max_vec, pch = 19, col = 3)



}

#plot phase dependent transport performance with connectivity trace
# or fixed value(Energy) trace

plot_performance_vs_trace = function(gr, by = .1)
{
  #convert to one element list if you get just a graph
  if(class(gr)== "graph")
    gr = list(gr)

	prog = seq(0,2*pi,by)

	performance_plot()

	codes = sapply(gr, "[[", "code")
	leg = c(paste(codes,"*e"), paste(codes,"*c"))
	legend(.2,1.3, legend = leg, col = 1:(2*length(gr)), pch = 1)

	gr = lapply(gr,retrace_E)
	max_vec = lapply(gr,performance,sample = prog)

	for(  i in seq_along(max_vec) )
    lines( prog, max_vec[[i]], col = i)

  gr = lapply(gr,retrace_conn)
	max_vec = lapply(gr,performance,sample = prog)

	for(  i in seq_along(max_vec) )
    lines( prog, max_vec[[i]], col = length(gr)+i)

}

# same as plot_performacne_vs_trace but the focus is arrival time
plot_t_performance_vs_trace = function(gr, by = .1)
{
  #convert to one element list if you get just a graph
  if(class(gr)== "graph")
    gr = list(gr)

	prog = seq(0,2*pi, by)
  t_max = max(sapply(gr,max_search_time))

	plot(-10,-10,ylim = c(0,t_max), xlim= c(0,2*pi),xlab = "e^ phi / (0 -> 2pi)", ylab="t_p_max")
  grid()

	codes = sapply(gr, "[[", "code")
	leg = c(paste(codes,"*e"), paste(codes,"*c"))
	legend(.2,t_max, legend = leg, col = 1:(2*length(gr)), pch = 1)

	gr = lapply(gr,retrace_E)
	max_vec = lapply(gr,t_performance,prog)

	for(  i in seq_along(max_vec) )
    lines( prog, max_vec[[i]], col = i)

  gr = lapply(gr,retrace_conn)
	max_vec = lapply(gr,t_performance, prog)

	for(  i in seq_along(max_vec) )
    lines( prog, max_vec[[i]], col = length(gr)+i)

}

# increasing size graph best performance with two trace generation method
plot_N_progression_vs_trace = function(code,range = 2:10)
{
  prog = seq(0,2*pi,.1)

  lg = c(paste(code,"N *e",sep = "") , paste(code,"N *c",sep = ""))

	plot(-10,-10,ylim = c(0,1.3), xlim= range(range),xlab = "N", ylab="p_max_phi")
	legend("topleft", legend = lg, col = 2:3, pch = 1)
  grid()

  gr_vec = get_N_topo(code ,range)

  ####
  gr_vec = lapply(gr_vec, retrace_E)
  opt_and_plot_list(gr_vec, mode = "sites", col = 2)

	####
  gr_vec = lapply(gr_vec, retrace_conn)
  opt_and_plot_list(gr_vec, mode = "sites", col = 3)

}

# increasing size graph best performance transport time with two trace generation method
plot_t_N_progression_vs_trace = function(code,range = 2:10)
{
  prog = seq(0,2*pi,.1)

  lg = c(paste(code,"N *e",sep = "") , paste(code,"N *c",sep = ""))
  t_max = max_search_time( ringNhandle(max(range)))

	plot(-10,-10,ylim = c(0,t_max), xlim= range(range),xlab = "N", ylab="t")
	legend("topleft", legend = lg, col = 2:3, pch = 1)
  grid()

  gr_vec = get_N_topo(code ,range)

  ####
  gr_vec = lapply(gr_vec, retrace_E)
  opt_and_t_plot_list(gr_vec, mode = "sites", col = 2)

	####
  gr_vec = lapply(gr_vec, retrace_conn)
  opt_and_t_plot_list(gr_vec, mode = "sites", col = 3)

}

#trasport rateo best P/ arrival time for increasing size 3 example topos
plot_N_transport_rateo = function(range = 2:10)
{
	prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,1), xlim= range(range),xlab = "N", ylab="otp_rateo")
	legend("topleft", legend = c("CN", "ChN","LN"), col = 1:3, pch = 1)
  grid()

	####
	gr_vec = lapply(range,ring)
	gr_vec = lapply(gr_vec,optimize_gr)

	gr_m = lapply(gr_vec, locate_max)
	maxima = sapply(gr_m, "[[", "m")
	times = sapply(gr_m, "[[", "t")
	rateos = maxima/times

	lines( range, rateos, pch = 19, col = 1)
	points( range, rateos, pch = 19, col = 1)

	####
	gr_vec = lapply(range,ringNhandle)
	gr_vec = lapply(gr_vec,optimize_gr)

	gr_m = lapply(gr_vec, locate_max)
	maxima = sapply(gr_m, "[[", "m")
	times = sapply(gr_m, "[[", "t")
	rateos = maxima/times

	lines( range, rateos, pch = 19, col = 2)
	points( range, rateos, pch = 19, col = 2)

	####
	gr_vec = lapply(range,line)
	gr_vec = lapply(gr_vec,optimize_gr)

	gr_m = lapply(gr_vec, locate_max)
	maxima = sapply(gr_m, "[[", "m")
	times = sapply(gr_m, "[[", "t")
	rateos = maxima/times

	lines( range, rateos, pch = 19, col = 3)
	points( range, rateos, pch = 19, col = 3)

}


#brute force check in the neighbouhr of a maximum
check_max = function( gr, t_tol = .001, phi_tol = .001, sub = 100)
{

  phi = optimum_phase(gr)
  phi_lm =optimum_phase_pm(gr)*pi

  #get scan grid around optimum phase
  phi_seq = seq( phi - phi_tol, phi + phi_tol, phi_tol/sub)
  phi_seq = exp(1i*pi*phi_seq)

  phi_lm_seq = seq( phi_lm - phi_tol, phi_lm + phi_tol, phi_tol/sub)
  phi_lm_seq = exp(1i*pi*phi_lm_seq)

  max = locate_max( optimize_gr(gr))
  max_lm = locate_max( optimize_gr(gr, mode = "lm"))

  #get scan grid around time of transport maximum
  t_seq = seq( max$t - t_tol, max$t + t_tol, t_tol/100)
  t_lm_seq = seq( max_lm$t - t_tol, max_lm$t + t_tol, t_tol/100)


  #brute force grid scan for local maximum

  fix_max = function(n_phi)
  {
    gr = rephase(gr,n_phi)

    p_vec = Mod(evo_l_t(gr,gr$start, t_seq))^2
    p_vec = p_vec[gr$target,]

    return ( c( max(p_vec), which.max(p_vec)) )
  }

  best_vec = sapply(phi_seq,fix_max)

  best_pos = which.max( best_vec[1,])
  best = best_vec[,best_pos]

  print( paste("Probabilità massima: ",best[1]))
  print( paste("T del massimo: ",t_seq[best[2]]))
  print( paste("fase_applicata : ",phi_seq[best_pos]))

  print(paste("Massimo algo p: " ,max$m," t:", max$t))
  print(paste("Fase algo", exp(1i*pi*phi)))
}

#optimize and plot on an already existing canvas the list probability progression
opt_and_plot_list = function(gr_list,opt_mode = "",mode = "range",range, nophi = FALSE,LINES = TRUE,POINTS = TRUE, ...)
{

  max_vec_nophi = sapply(gr_list,p_max_func)

  gr_list = lapply(gr_list,optimize_gr, mode = opt_mode)

  max_vec = sapply(gr_list,p_max_func)

  if(mode != "range")
    x = get_list_x(gr_list, mode)
  else
    x = range


  if(LINES)
    lines( x, max_vec, ...)
  if(nophi)
    lines( x, max_vec_nophi, lty = "dotted", ...)
  if(POINTS)
    points( x, max_vec, ...)

	return( max(max_vec))
}

#optimize and plot on an already existing canvas the list time of arrival progression
opt_and_t_plot_list = function(gr_list,opt_mode = "",mode = "range",range, nophi = FALSE,LINES = TRUE, POINTS = TRUE, ...)
{
  max_vec_nophi = sapply(gr_list,t_p_max_func)

  gr_list = lapply(gr_list,optimize_gr, mode = opt_mode)

  max_vec = sapply(gr_list,t_p_max_func)

  if(mode != "range")
    x = get_list_x(gr_list, mode)
  else
    x = range

  print(length(x))
  print(length(gr_list))

  if(LINES)
    lines( x, max_vec, ...)
  if(nophi)
    lines( x, max_vec_nophi, lty = "dotted", ...)
  if(POINTS)
    points( x, max_vec, ...)

	return("time")
}

#reference utility line for performance trends
plot_line_trend = function(range, mode = "dist", col = 3)
{
  gr_list = lapply(range,line)

  opt_and_plot_list(gr_list, mode = mode, opt_mode = "none",POINTS = F,col = col)
}

performance_plot = function( title = "")
{
  plot(-10,-10,
       ylim = c(0,1.3), xlim= c(0,2*pi),
       xlab = expression(phi), ylab=expression("p"["max"]))
  title(main =title)

  grid()
}

t_performance_plot = function( title = "", t_max)
{
  plot(-10,-10,
       ylim = c(0,t_max), xlim= c(0,2*pi),
       xlab = expression(phi), ylab=expression(t[p["max"]]))
  title(main =title)

  grid()
}

progression_plot = function( title = "", xmax = 20,x_lab = "N")
{

  plot(-10,-10,
       ylim = c(0,1.3), xlim= c(0,xmax),
       xlab = x_lab, ylab=expression("p"["max"]))
  title(main =title)

  grid()
}

t_progression_plot = function( title = "",ymax=10, xmax = 20,x_lab = "N")
{

  plot(-10,-10,
       ylim = c(0,ymax), xlim= c(0,xmax),
       xlab = x_lab, ylab=expression("t"["p"["max"]]))
  title(main =title)

  grid()
}



