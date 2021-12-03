# decompose l th localized state into matrix eigenvectors

decompose_localized = function(l)
{
	H_ev = eigen(H)

	ret = Conj(H_ev$vectors[l,])

	ret	
}

# recompose eigenvector basis amplitudes into 
# localized states amplitude

recompose = function(eigen_A)
{
	H_ev$vectors %*% eigen_A
}


evo_phase = function(t)
{
	ret = function(dt) exp(-1i*H_ev$values*dt)
	drop(sapply(t,ret))
	
}

evo_l_t = function(l, t)
{
	A_t = decompose_localized(l)*evo_phase(t)
	A_t = recompose(A_t)

	A_t
}

evo_full = function(init = 1,from = 0, to=5, by = .1 )
{
	t_seq = seq(from, to, by)
	evo_l_t(l=init,t_seq)
}

# plot time evulution with different color 
# for different sites for an arbitrary complex matrtx

plot_evo_mat = function(init=1,from = 0, to=5, by = .1)
{
	prepare_globals()	
	size = dim(H)[1]
	p_vec = Mod(evo_full(init, from, to ,by))^2

	plot(0,0,ylim = c(0,1), xlim= c(from,to),xlab = "t", ylab="p(t)")
	legend(3,1, legend = 1:size , col = 1:size   ,lty = 1)

	forplo(a in 1:size )
	lines(seq(from,to,by),p_vec[a,]/sum(p_vec[,1]), col= a, pch=19)


	#debug stuff, works only with a even N ring
	max(p_vec[size,])
	
}

# plot the probability of transport between 
# leftmost site and target(rightmost) site
# in a N sites chain and  N sites + 2 handles

plot_evo_vs_chain = function(N, phi = 0,from = 0 ,to= 10,by = .1)
{	
	
	plot(0,0,ylim = c(0,1), xlim= c(from,to),xlab = "t", ylab="p(t)")
	legend(3,1, legend = c(gsub("N", N, "C(N)"),gsub("N", N, "Ch(N)") ), col = 2:3   ,lty = 1)

	assign("m", adjacency(N), envir = .GlobalEnv)
	assign("m", rephase(phi), envir = .GlobalEnv)
	prepare_globals()
	target = N%/%2+1

	p_vec = Mod(evo_full(init = 1, from, to ,by))^2
	paste("Max loop:" , max(p_vec[target,]) )

	lines(seq(from,to,by),p_vec[target,]/sum(p_vec[,1]), col= 2, pch=19)


	assign("m", chainNhandle(N), envir = .GlobalEnv)
	assign("m", rephase(phi), envir = .GlobalEnv)
	prepare_globals()
	target = N+2

	p_vec = Mod(evo_full(init = N+1, from, to ,by))^2
	paste("Max chain:" , max(p_vec[target,]) )

	lines(seq(from,to,by),p_vec[target,]/sum(p_vec[,1]), col= 3, pch=19)


}

plot_dist_mat = function(init=1,from = 0, to=5, by = 1)
{
	t_seq = seq(from,to,by)
	prepare_globals()
	size = dim(H)[1]
	p_vec = Mod(evo_full(init, from, to ,by))^2

	plot(0,0,ylim = c(0,1), xlim= c(1,size ),xlab = "l", ylab="p(l)")

	for(a in 1:length(t_seq))
	  lines(1:size, p_vec[,a]/sum(p_vec[,1]), col= a, pch=19)

	legend(3,1, legend = t_seq, col = 1:length(t_seq)  ,lty = 1)
	
}