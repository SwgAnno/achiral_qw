#all the support methods assume the existence of a variable gr
# which stores a Graph object(see objMat.R)

# decompose l th localized state into matrix eigenvectors

decompose_localized = function(gr,l)
	return(  Conj(gr$eigen$vectors[l,]))

# recompose eigenvector basis amplitudes into
# localized states amplitude

recompose = function(gr,eigen_A)
	return(  gr$eigen$vectors %*% eigen_A )


evo_phase = function(gr,t)
{
	exp_map = function(dt) exp(-1i* gr$eigen$values *dt)
	return(  drop(sapply(t,exp_map)))
}


evo_l_t = function(gr, l = gr$start, t)
{
	A_t = decompose_localized(gr,l)*evo_phase(gr,t)
	return(  recompose(gr,A_t))
}


evo_l_lm = function(gr,l = gr$start, model = get_topo_lm(gr))
{
  t = get_lm_time(gr,model)

  return(evo_l_t(gr, l, t))
}

evo_full = function(gr, from = 0, to=5, by = .1 )
{
	t_seq = seq(from, to, by)
	return(  evo_l_t(gr, l= gr$start,t_seq))
}

evo_full_target = function(gr, target = gr$target, from = 0, to = max_search_time(gr), by = .1 )
{
	t_seq = seq(from, to, by)
	p_mat = evo_l_t(gr, l= gr$start,t_seq)

	return(p_mat[target,])
}

#derivative of site probability l at time t
deriv_p_l_t = function(gr,t, l = gr$target)
{
  A_t = decompose_localized(gr,gr$start)*evo_phase(gr,t)

  ret = Conj(gr$eigen$vectors) %*%(gr$eigen$values*Conj(A_t))
  ret = ret* recompose(gr,A_t)

  return(  -2*Im(ret[l,]))

}

deriv_full = function(gr, from = 0, to=5, by = .1 )
{
	t_seq = seq(from, to, by)
	return(  deriv_p_l_t(gr,t_seq))
}



# plot time evulution with different color
# for different sites for an arbitrary complex matrtx

plot_evo_mat = function(gr,from = 0, to= max_search_time(gr), by = .1)
{

	size = gr$N
	p_vec = Mod(evo_full(gr,from, to ,by))^2

	plot(0,0,ylim = c(0,1), xlim= c(from,to),xlab = "t", ylab="p(t)")
	legend(3,1, legend = 1:size , col = 1:size   ,lty = 1)
  grid()

	for(a in 1:size )
	lines(seq(from,to,by),p_vec[a,]/sum(p_vec[,1]), col= a, pch=19)


	#debug stuff, tells $target site maximum prob
	 print(max(p_vec[gr$target,]))

}

#restrict evo_mat to display just an handful of target sites
# initialized wiht just the target one
plot_evo_mat_target = function(gr,sites = gr$target,from = 0, to= max_search_time(gr), by = .1)
{

	size = gr$N
	p_vec = Mod(evo_full(gr,from, to ,by))^2

	plot(0,0,ylim = c(0,1), xlim= c(from,to),xlab = "t", ylab="p(t)")
	legend(3,1, legend = sites , col = sites   ,lty = 1)
  grid()

	for(a in seq_along(sites) )
	lines(seq(from,to,by),p_vec[sites[a],]/sum(p_vec[,1]), col= sites[a], pch=19)


	#debug stuff, tells $target site maximum prob
	 print(max(p_vec[gr$target,]))

}

#see how target site probability evolution changes with phase

#with higher number of phases uses diagonal rephase
plot_evo_mat_target_phase = function(gr, from = 0, to= max_search_time(gr), by = .1, phase_by = .1, SURFACE = FALSE)
{
  prog = seq(0,2*pi, phase_by)
  n = length(gr$re_coord)
  t_seq = seq(from,to,by)

  phi_mat = matrix(  rep(exp(1i*prog), n), ncol = n)

  gr_vec = apply(phi_mat,1, rephase, gr = gr)

  p_mat = sapply(gr_vec, evo_full_target, target = gr$target, from, to, by)

  #print(p_mat)
  p_mat = Mod(p_mat)^2

  colnames(p_mat) = prog
  rownames(p_mat) = t_seq

  if(SURFACE)
    rotate.persp( t_seq, prog, p_mat)
  else
    image.plot(t_seq, prog, p_mat, xlab = "t", ylab = expression(phi))


}

#plot p value of target site for a ring of N site with and without one site handles
# graphs are rephased once with the same phase
plot_evo_vs_ring = function(N, phi = 0,from = 0 ,to= max_search_time(ring(N)),by = .1)
{

	plot(0,0,ylim = c(0,1), xlim= c(from,to),xlab = "t", ylab="p(t)")
	legend(3,1, legend = nsub(N, c("C(N)","Ch(N)")) , col = 2:3   ,lty = 1)
  grid()


	gr = ring(N)
	gr = rephase(gr,phi)

	p_vec = Mod(evo_full(gr, from, to ,by))^2
	print(paste("Max ring:" , max(p_vec[gr$target,]) ))

	lines(seq(from,to,by),p_vec[gr$target,]/sum(p_vec[,1]), col= 2, pch=19)


	gr = ringNhandle(N)
	gr = rephase(gr,phi)

	p_vec = Mod(evo_full(gr, from, to ,by))^2
	print (paste("Max ring+handles:" , max(p_vec[gr$target,]) ))

	lines(seq(from,to,by),p_vec[gr$target,]/sum(p_vec[,1]), col= 3, pch=19)


}

#same as plot_evo_vs_ring but graphs are shown with their own best rephase
plot_fair_evo_vs_ring = function(N, phi = 0,from = 0 ,to= max_search_time(ring(N)),by = .1)
{

	plot(0,0,ylim = c(0,1), xlim= c(from,to),xlab = "t", ylab="p(t)")
	legend(3,1, legend = nsub(N, c("C(N)","Ch(N)")), col = 2:3   ,lty = 1)
  grid()

	gr = ring(N)
	gr = optimize_gr(gr)

	p_vec = Mod(evo_full(gr, from, to ,by))^2
	print(paste("Max ring:" , max(p_vec[gr$target,]) ))

	lines(seq(from,to,by),p_vec[gr$target,]/sum(p_vec[,1]), col= 2, pch=19)


	gr = ringNhandle(N)
	gr = optimize_gr(gr)

	p_vec = Mod(evo_full(gr, from, to ,by))^2
	print(paste("Max ring+handles:" , max(p_vec[gr$target,]) ))

	lines(seq(from,to,by),p_vec[gr$target,]/sum(p_vec[,1]), col= 3, pch=19)


}

# plot the distubution of probability sampled at some imtervals
# protip: it looks awful
plot_dist_mat = function(gr,from = 0, to=5, by = 1)
{
	t_seq = seq(from,to,by)
	size = gr$N
	p_vec = Mod(evo_full(gr, from, to ,by))^2

	plot(0,0,ylim = c(0,1), xlim= c(1,size ),xlab = "l", ylab="p(l)")
  grid()

	for(a in 1:length(t_seq))
	  lines(1:size, p_vec[,a]/sum(p_vec[,1]), col= a, pch=19)

	legend(3,1, legend = t_seq, col = 1:length(t_seq)  ,lty = 1)

}

# plot the probability evolution of every site in a heatmap fashion
plot_heat_mat = function(gr,from = 0, to = max_search_time(gr), by = .1)
{
	t_seq = seq(from,to,by)
	size = gr$N
	p_vec = t (Mod(evo_full(gr, from, to ,by))^2 )

	colnames(p_vec) = 1:gr$N
	rownames(p_vec) = t_seq

  plot(0,0,ylim = c(0,1), xlim= c(1,size ),xlab = "l", ylab="p(l)")
  image.plot(t_seq, 1:gr$N, p_vec, xlab = "t", ylab= "l")
}

# plot the probability evolution of every site in a surface mesh fashion
plot_surf_mat = function(gr,from = 0, to= max_search_time(gr), by = .1)
{
	t_seq = seq(from,to,by)
	size = gr$N
	p_vec = t (Mod(evo_full(gr, from, to ,by))^2 )

	colnames(p_vec) = 1:gr$N
	rownames(p_vec) = t_seq

  plot(0,0,ylim = c(0,1), xlim= c(1,size ),xlab = "l", ylab="p(l)")
  #rotate.persp(t_seq, 1:gr$N, p_vec, xlab = "t", ylab= "l")
  persp3D(t_seq, 1:gr$N, p_vec,
          axes=TRUE,box=TRUE,
          nticks=5,ticktype="detailed",
          xlab = "t", ylab = "l", zlab = "p")
}

