#minmax

deriv_p_l_t = function(init = 1,l = 1,t)
{
  A_t = decompose_localized(init)*evo_phase(t)

  ret = Conj(H_ev$vectors) %*%(H_ev$values*Conj(A_t))
  ret = ret* recompose(A_t)

  -2*Im(ret[l,])

}

deriv_roots = function(init = 1,l=1, from = 0,to =5)
{

	f = function(t) deriv_p_l_t(init,l,t)
	uniroot.all(f,lower = from, upper = to)
	
}

p_max = function(init = 1,l=1, from = 0,to =10, by = .1)
{
	candidates = deriv_roots(init,l,from,to)
	p_vec = Mod(evo_l_t(init, candidates))^2

	return( max(p_vec[l,]))
}

evo_vs_deriv = function(init=1,l=1,from = 0, to=5, by = .1)
{
	t_seq = seq(from,to,by)
	prepare_globals()

	p_vec = Mod(evo_full(init, from, to ,by))^2
	deriv_vec = deriv_p_l_t(init=init, l=l,t = t_seq)

	plot(0,0,ylim = c(-2,2), xlim= c(from,to ),xlab = "l", ylab="p(l)")

	lines(t_seq, p_vec[1,]/sum(p_vec[,1]), col= l)
	lines(t_seq, deriv_vec,col = l+1)

	legend(3,1, legend = c("P(t)","P'(t)"), col = c(l,l+1)  ,lty = 1)
}

performance = function(sample, init = 1, l)
{
	phi_vec = exp(sample *1i)
	ret = c(length(sample))

	for(i in 1:length(phi_vec))
	{
		assign("m" ,rephase( phi_vec[i]),envir = .GlobalEnv)

		prepare_globals()
		ret[i] = p_max(init = init ,l=l)
	}

	return(ret)
}
	
plot_performance = function()
{
	prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,1.3), xlim= c(0,2*pi),xlab = "e^ phi / (0 -> 2pi)", ylab="p_max")
	
	labels = c("C(7)","Ch(7)")
	legend(.2,1.3, legend = labels, col = 2:3, pch = 1)

	assign("m" , adjacency(7),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog, l = N%/%2+1)

	lines( prog, max_vec, col = 2)

	assign("m" , chainNhandle(7),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog,init = N+1 ,l = N+2)

	lines( prog, max_vec, col = 3)

}

plot_even_progression = function()
{
	prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,1.3), xlim= c(0,2*pi),xlab = "e^ phi / (0 -> 2pi)", ylab="p_max")
	
	labels = c("C(4)","C(6)", "C(8)","C(10)")
	legend(.2,1.3, legend = labels, col = 2:5, pch = 1)

	assign("m" , adjacency(4),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog, l = N%/%2+1)

	lines( prog, max_vec, col = N/2)


	assign("m" , adjacency(6),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog, l = N%/%2+1)

	lines( prog, max_vec, col = N/2)


	assign("m" , adjacency(8),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog, l = N%/%2+1)

	lines( prog, max_vec, col = N/2)


	assign("m" , adjacency(10),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog, l = N%/%2+1)

	lines( prog, max_vec, col = N/2)



}

plot_vs_chain = function(N)
{
	prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,1.3), xlim= c(0,2*pi),xlab = "e^ phi / (0 -> 2pi)", ylab="p_max")
	
	labels = c(gsub("N", N, "C(N)"),gsub("N", N, "Ch(N)"))
	legend(.2,1.3, legend = labels, col = 2:3, pch = 1)

	assign("m" , adjacency(N),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog, l = N%/%2+1)

	lines( prog, max_vec, col = 2)

	assign("m" , chainNhandle(N),envir = .GlobalEnv)
	prepare_globals()

	max_vec = performance(sample = prog,init = N+1 ,l = N+2)

	lines( prog, max_vec, col = 3)

}

plot_N_progression = function(range = 2:10)
{
	prog = seq(0,2*pi,.1)

	plot(-10,-10,ylim = c(0,1.3), xlim= range[c(1,length(range))],xlab = "N", ylab="p_max_phi")
	legend(6,1.2, legend = c("CN", "ChN"), col = 1:2, pch = 1)

	max_vec = c(length(range))

	for( i in 1: length(range))
	{
		N = range[i]
		assign("m" , adjacency(N),envir = .GlobalEnv)
		prepare_globals()

		max_N_vec = performance(sample = prog, l = N%/%2+1)

		max_vec[i] = max( max_N_vec)

	}

	lines( range, max_vec, pch = 19, col = 1)
	points( range, max_vec, pch = 19, col = 1)

	for( i in 1: length(range))
	{
		N = range[i]
		assign("m" , chainNhandle(N),envir = .GlobalEnv)
		prepare_globals()

		max_N_vec = performance(sample = prog, init = N+1,l = N+2)

		max_vec[i] = max( max_N_vec)

	}

	lines( range, max_vec, pch = 19, col = 2)
	points( range, max_vec, pch = 19, col = 2)



}