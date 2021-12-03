#you have to define N (lenght of the ring) in order to
#have all those function work


#bloch wave number and energy value 

bloch_k = function(n)
{
	pi*2*n/N
}

bloch_E = function(n)
{
	2-2*cos(bloch_k(n))
}

#probability of being in site l at time t starting from 0 at time 0

p_l_t = function(pos,t)
{

	i = 0:(N-1)
	result =  exp( (bloch_k(i)*pos - bloch_E(i)*t)*1i)

	Mod(sum(result))^2
}

evolution_l = function(l, from = 0, to=1, by = .1)
{
	t = seq(from,to, by)
	p_vec = sapply(t, p_l_t, pos=l)

	p_vec
}

state_t = function(time)
{
	l = (-N%/%2):(N%/%2)
	p_vec = sapply(l, p_l_t, time)

	p_vec
}

#plot time evulution with different color for different sites

plot_evo = function(from = 0, to=5, by = .1)
{

	plot(0,0,ylim = c(0,1), xlim= c(from,to),xlab = "t", ylab="p(t)")

	i = 0:(N-1)
	p_vec = sapply(i,evolution_l,from,to,by)
	
	for(a in 1:N)
	lines(seq(from,to,by),p_vec[,a]/sum(p_vec[1,]), col= a, pch=19)
	
}

# plot the site probability distribution with different
# color for increasing time

plot_dist = function(from = 0, to=5, by = 1)
{

	plot(0,0,ylim = c(0,1), xlim= c(-N%/%2,N%/%2),xlab = "l", ylab="p(l)")

	i = seq(from,to,by)
	p_vec = sapply(i,state_t)
	
	for(a in 1:((to-from)/by))
	lines( (-N%/%2):(N%/%2) ,p_vec[,a]/sum(p_vec[,1]), col= a, pch=19)
	
}

plot_dist_grey = function(from = 0, to=5, by = 1)
{

	plot(100,100,ylim = c(from,to), xlim= c(-N%/%2,N%/%2),xlab = "l", ylab="p(l)")

	i = seq(from,to,by)
	size = (to-from)/by
	p_vec = sapply(i,state_t)
	p_vec = p_vec/sum(p_vec[,1]) #normalize

	palette = colorRamp(c("white","red"))
	prob_color = palette(seq(0,1,.01))
	prob_color = prob_color[ round(p_vec*100)+1]	
	dim(prob_color) = dim(p_vec)
	
	for(a in 1:((to-from)/by))
	points( (-N%/%2):(N%/%2) ,rep(i[a],N+1), col= prob_color[,a], pch=19)
	
}