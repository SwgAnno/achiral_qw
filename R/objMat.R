#matrix object


# class mat()	// it's a list atm

# int N 			// number of sites
# char* code  // plot code (utility)
# comp[N][N] mat 		// adjaacency matrix
# int[X][2] re_coord	// rephase coord
#				   ( just one of the two symetric matrix entry)
# eigen			// eigenvalue|vectors list
# int start			// start site for evolution
# int target 		// target site ''

# "constructor" for a N sites|no connection graph

# 0 size graph are defined but usually ill formed

plot.graph = function(gr,...)
{
  plot( network_from_graph(gr),...)
}

graph = function(n_N = 4)
{
	out = list()
	out$N = n_N
	out$code = "e"

	out$mat = matrix(complex(n_N*n_N), ncol = n_N)

	if(n_N == 0)
    return(out)

	for(i in 1:n_N)
		out$mat[i,i] = 1

	out$eigen = eigen(out$mat)

	out$re_coord = list(c(-1,-1))
	out$start = 1
	out$target = 1

	class(out) <- "graph"

	return(out)

}

# returns N sites ring
ring = function(N, E=2,conn_tr = FALSE)
{
	out = graph(N)
	out$code = paste("C",N,sep="")

	if(N == 0)
    return(out)

	for( i in 1:N)
	{
 		out$mat[i,(i%%N+1)] <- as.complex(-1)
		out$mat[(i%%N+1),i] <- as.complex(-1)
  }

  if(conn_tr)
    out = retrace_conn(out)
  else
    out = retrace_E(out,E)

	out$eigen = eigen(out$mat)

	out$start = 1
	out$target = N%/%2+1

	if(N!=1)
    out$re_coord[[1]] <- c(1,N)


	return(out)
}

ringNhandle = function(N, E=2,conn_tr = FALSE)
{
	out = graph(N+2)
  out$code = paste("Ch",N,sep="")

  if(N == 0)
    return(out)

	for( i in 1:N)
	{
		out$mat[i,(i%%N+1)] <- as.complex(-1)
		out$mat[(i%%N+1),i] <- as.complex(-1)
  }

	out$mat[N+1,N+1] = 2
	out$mat[N+2,N+2] = 2

	out$mat[N+1,1] = -1
	out$mat[1,N+1] = -1

	out$mat[N+2,N%/%2+1] = -1
	out$mat[N%/%2+1,N+2] = -1

  if(conn_tr)
    out = retrace_conn(out)
  else
    out = retrace_E(out,E)

	out$eigen = eigen(out$mat)

	out$start = N+1
	out$target = N+2

	if(N!=1)
    out$re_coord[[1]] <- c(1,N)

	return(out)
}

line = function(N, E=2,conn_tr = FALSE)
{
	out = graph(N)
  out$code = paste("L",N,sep="")

  if(N == 0)
    return(out)

	for( i in 1:N-1)
	{
 		out$mat[i,(i+1)] <- as.complex(-1)
		out$mat[(i+1),i] <- as.complex(-1)
  }

  if(conn_tr)
    out = retrace_conn(out)
  else
    out = retrace_E(out,E)

	out$eigen = eigen(out$mat)

	out$start = 1
	out$target = N
	out$re_coord[[1]] <- c(-1,-1)


	return(out)
}


## half cut ===> cut = c(1+ (N+3)%/%4, N- (N-3)%/%4

ringNcut = function(N,cut = c(1,1+ N%/%2 ), HANDLES = TRUE, E=2,conn_tr = FALSE)
{


  if(HANDLES)
    out = ringNhandle(N, E, conn_tr)
  else
    out = ring(N, E, conn_tr)

  if(N<=3)
    return(out)

  out$code = nsub(N,"RcN")

  if(cut[1]!=cut[2])
  {
    #add the cut
    out$mat[cut[1],cut[2]] <- as.complex(-1)
    out$mat[cut[2],cut[1]] <- as.complex(-1)

    #reposition re_coord links at the ends of the ring
    out$re_coord[[1]] <- c(1,2)
    out$re_coord[[2]] <- c(cut[1],cut[2])
  }

  out$eigen = eigen(out$mat)

  return(out)

}

#create a chain with gr as a unit
chain = function(gr, rep, space = 0, HANDLES = TRUE)
{
  if(HANDLES)
    base_gr = join_link( line(1), gr)
  else
    base_gr = gr

  if(rep != 1)
  {
    unit = join_nolink(line(space+1), gr)

    for(i in 1:(rep-1))
      base_gr = join_nolink(base_gr, unit)

  }

  if(HANDLES)
    base_gr = join_link( base_gr, line(1))

  return(base_gr)

}

reverse = function(gr)
{
  s = gr$start
  t = gr$target

  gr$start = t
  gr$target = s

  return(gr)
}

topo_from_code = function( code)
  return (str_extract(code, "[a-zA-Z]+")  )

size_from_code = function(code)
  return (as.numeric(str_extract(code, "\\d+"))  )


graph_from_code = function( codes )
{
  # codes are supposed to be "ccd+"
  # two char identifier and max two digit size

  # "Cd"    --> ring with d sites
  # "Chd"   --> ring with d sites + 2 handles
  # "Ld"    --> Line with d sites

  out = list()

  for( i in seq_along(codes))
  {
    cur = codes[i]

    gr_code = str_extract_all(cur, "[a-zA-Z]+\\d+")[[1]]
    op_code = str_extract_all(cur, "[+/]")[[1]]

    gr_vec = graph_from_mono_code(gr_code)

    tot_gr = gr_vec[[1]]

    if(str_detect(cur, "[+/]"))
    {
      for(m in seq_along(op_code) )
      {
        if(op_code[m] == "+")
          tot_gr = join_link(tot_gr, gr_vec[[m+1]])

        else if(op_code[m] == "/")
          tot_gr = join_nolink(tot_gr,gr_vec[[m+1]])

      }
    }

    out[[i]]= tot_gr
  }

  return(out)
}

gfc = function(codes) return(graph_from_code(codes)[[1]])

graph_from_mono_code = function(codes)
{
  # codes are supposed to be "ccd+"
  # two char identifier and max two digit size

  # "Cd"    --> ring with d sites
  # "Chd"   --> ring with d sites + 2 handles
  # "Ld"    --> Line with d sites
  # "Rcd"   --> Ring with d sites and a cut

  out = list()

  sizes = sapply(codes,size_from_code)
  topos = sapply(codes,topo_from_code)
  for(i in seq_along(codes))
  {
    if(topos[i] == "C")
      out[[i]] = ring(sizes[i])

    else if( topos[i] == "Ch")
      out[[i]] = ringNhandle(sizes[i])

    else if( topos[i] == "L")
      out[[i]] = line(sizes[i])

    else if( topos[i] == "Rc")
      out[[i]] = ringNcut(sizes[i])

    else
    {
      print("Error in topo from mono code, topology not found")
      out = line(0)
    }

  }

  return(out)
}

max_search_time = function(gr)
{
  return( gr$N * getTC())
}


#rephase connections listed in graph object
rephase = function(gr , phi = 1i)
{

  #assert_is_complex(phi)

  if( length(gr$re_coord) != length(phi))
  {
    print("rephase() error, wrong number of given phases")
    return(gr)
  }

  if( gr$re_coord[[1]][1] == -1)
    return(gr)


	for( i in 1:length(gr$re_coord))
	{
		v = gr$re_coord[[i]]
		gr$mat[v[1],v[2]] = -1*phi[i]
		gr$mat[v[2],v[1]] = -1*Conj(phi[i])


	}

	gr$eigen = eigen(gr$mat)

	return(gr)
}

#redefines trace looking for connectivity
retrace_conn = function(gr)
{
  for(i in 1:gr$N)
  {
    count = 0
    for(l in 1:gr$N)
    {
      if(gr$mat[i,l]!=0 && i!=l)
        count = count+1

    }

    gr$mat[i,i] = count
  }


  gr$eigen = eigen(gr$mat)

  return(gr)
}


#redefines trace with a default value of energy
retrace_E = function(gr,E = 2)
{
  for(i in 1:gr$N)
    gr$mat[i,i] = E

  gr$eigen = eigen(gr$mat)

  return(gr)

}

#join two graph with a link between target1 e start2
join_link = function(gr1, gr2)
{
  if( gr1$N == 0)
    return(gr2)
  else if( gr2$N ==0)
    return(gr1)

  oN = gr1$N + gr2$N
  out = graph( oN )

  out$code = paste(gr1$code, "+", gr2$code, sep = "")

  #copy old mats
  out$mat[1:gr1$N,1:gr1$N] = gr1$mat
  out$mat[(gr1$N+1):oN, (gr1$N+1):oN] = gr2$mat

  #add new link
  out$mat[gr1$target, gr1$N+ gr2$start] = as.complex(-1)
  out$mat[gr1$N+ gr2$start, gr1$target] = as.complex(-1)

  out$eigen = eigen(out$mat)

	out$start = gr1$start
	out$target = gr1$N +gr2$target

  #shift re_coord of second graph
  rec2 = list()

  for(i in seq_along(gr2$re_coord))
  {
    if(gr2$re_coord[[i]][1] == -1 )
      rec2[[i]] = gr2$re_coord[[i]]
    else
      rec2[[i]] = gr2$re_coord[[i]]+ gr1$N
  }


	out$re_coord = c( gr1$re_coord, rec2)

  return( clean_re_coord(out) )
}

#join two graph, now target1 and start2 are the same site
join_nolink = function(gr1, gr2)
{
  if( gr1$N == 0 || gr1$N == 1)
    return(gr2)
  else if( gr2$N == 0 || gr2$N == 1)
    return(gr1)

  oN = gr1$N + gr2$N -1
  out = graph( oN )

  out$code = paste(gr1$code, "/", gr2$code, sep = "")

  #copy old mat1
  out$mat[1:gr1$N,1:gr1$N] = gr1$mat

  # copy mat2 without the site that's gonna be joined
  out$mat[(gr1$N+1):oN, (gr1$N+1):oN] = gr2$mat[-gr2$start,-gr2$start]

  # add to-be-joined site links
  out$mat[gr1$target, (gr1$N+1):oN] = gr2$mat[gr2$start,-gr2$start]
  out$mat[(gr1$N+1):oN, gr1$target] = gr2$mat[-gr2$start,gr2$start]

  out$eigen = eigen(out$mat)

	out$start = gr1$start
	out$target = gr1$N +gr2$target-1

  #shift re_coord of second graph

  #lookup array with changed site numbers
	new_n = c( 1:gr2$start, gr2$start:(gr2$N-1) )
  new_n = new_n + gr1$N
	new_n[gr2$start] = gr1$target

	#print(new_n)

	rec2 = list()

  for(i in seq_along(gr2$re_coord))
  {
    if(gr2$re_coord[[i]][1] == -1 )
      rec2[[i]] = gr2$re_coord[[i]]
    else
      rec2[[i]] = new_n[ gr2$re_coord[[i]] ]
  }

	out$re_coord = c( gr1$re_coord, rec2)

  return( clean_re_coord(out) )
}

add_handles = function( base, l = 1, mode = "both", fix = 0)
{
  if(mode == "both")
  {
    dx = l
    sx = l
  }
  if(mode == "sx")
  {
    dx = 0
    sx = l
  }
  if(mode == "dx")
  {
    dx = l
    sx = 0
  }
  if(mode == "fixl")
  {
    dx = l
    sx = fix
  }
  if(mode == "fixr")
  {
    dx = fix
    sx = l
  }
  out = join_link(line(sx),base)
  out = join_link(out,line(dx))

  return(out)
}

clean_re_coord = function(gr)
{

  toDel = function( a )
    return ( a[1] == -1)

  del = sapply( gr$re_coord, toDel)
  gr$re_coord = gr$re_coord[!del]

  if( length(gr$re_coord) == 0)
    gr$re_coord[[1]] = c(-1,-1)


  return ( gr )
}

# get best transport-wise phase according to optimization mode

optimize_gr = function(gr, mode = "", model = get_topo_lm(gr), by = .1)
{
  #print(paste("optimizing ", gr$code))

  if(mode == "")
    return(  rephase(gr, exp(1i*pi* optimum_phase(gr)) )  )

  else if(mode == "yolo")
    return(  rephase(gr, exp(1i*pi* optimum_phase_yolo(gr, by)) )  )

  else if(mode == "lm")
    return(  rephase(gr, exp(1i*pi* optimum_phase_lm(gr, model)) )  )

  else if(mode == "diag")
  {
    phi = optimum_phase_diag(gr)
    return(  rephase(gr, exp(1i*pi* rep(phi, length(gr$re_coord)) )  ))
  }
  else if(mode == "smart")
    return(  rephase(gr, exp(1i*pi* optimum_phase_smart(gr)) )  )

  else
    return(gr)


}

#n dof graph optimization
optimum_phase = function(gr)
{
  n = length( gr$re_coord)

  if(n==1)
    return(optimum_phase_1(gr))

  get_perf = function(phi)
    p_max( rephase(gr, exp(1i*phi)  ))

  if(GLOBAL_FIRST_MAX)
  get_perf = function(phi)
    p_max_first( rephase(gr, exp(1i*phi)  ))

  #CAUTION : number of start point shall scale with ^n but witn 10 dimentions is pure madness
  init = matrix( runif(GLOBAL_RANDOM_START_COUNT*n) ,ncol = n)
  init = init * 2*pi

  #optim control parameters
  my_control = list()
    my_control$fnscale = -1


  max = apply(init,1 ,optim, get_perf, method = "L-BFGS-B", lower = 0,upper = 2*pi, control = my_control)

  max_phi = sapply( max, "[[", "value")
  max_i = which.max(max_phi)

  return(max[[max_i]]$par /pi )

}

#1 dof graph optimization
#displays maximum in pi unit
#as search for optimum is done between 0 and 2 pi
#WARING : optimize can do strange stuff,double-check results
optimum_phase_1 = function(gr)
{
  get_perf = function(phi)
    p_max_func( rephase(gr, exp(1i*phi)  ))

  result = optimize(get_perf, interval = c(-2*pi,2*pi),maximum = TRUE,tol=0.0001)

  return( result$maximum/pi  )

}

optimum_phase_lm = function(gr, model = get_topo_lm(gr))
{
  get_perf = function(phi)
    p_lm( rephase(gr,exp(1i*phi)), model)

  result = optimize(get_perf, interval = c(-2*pi,2*pi),maximum = TRUE,tol=0.0001)

  return( result$maximum/pi  )
}

optimum_phase_yolo = function(gr, by = .1)
{
  prog = seq(0,2*pi,.1)
  perf = performance(sample = prog,gr)

  pos_max = which.max(perf)

  return( prog[pos_max]/pi )
}

optimum_phase_diag = function(gr)
{
  N_phi = length(gr$re_coord)

  get_perf = function(phi)
    p_max_func( rephase( gr, exp(1i*rep(phi,N_phi)) ) )

  result = optimize(get_perf, interval = c(-2*pi,2*pi),maximum = TRUE,tol=0.0001)

  return( result$maximum/pi  )
}

#find the best performance on the grid of frequent optimal values, that is the multiple of pi/2
# the size of p_mad grows as 4^N_phi, so it's a bit impractical for longer chains
optimum_phase_smart = function(gr)
{
  good_phases = c(0, pi/2, pi, 3/2*pi)
  N_phi = length(gr$re_coord)

  phase_list = rep( list(good_phases), N_phi)
  phase_grid = expand.grid(phase_list)

  #print(phase_grid)

  get_perf = function(phi)
    p_max_func( rephase( gr, exp(1i*phi) ) )

  p_vec = apply(phase_grid, 1,get_perf)
  result = as.numeric(phase_grid[which.max(p_vec),]/pi)

  #assert_is_numeric(result)
  return( result )
}

#return array of code topo graph with target range topology
get_N_topo = function(code,range)
{

	if(code=="C")
    return( lapply(range,ring))
  else if( code == "Ch")
    return( lapply(range,ringNhandle))
  else if(code == "L")
    return( lapply(range,line))

}

#optimize get_N_topo_result
# NOT GENERAL, only works for 1D graph
get_N_topo_opt = function(code, range)
{
  out = get_N_topo(code,range)

  if(code != "L")
  out = lapply(out,optimize_gr)

  return(out)
}


get_topo_lm = function(gr)
{
  topo = topo_from_code(gr$code)

  if(topo == "C")
    return( ring_t_lm())
  if(topo == "Ch")
    return( ringNh_t_lm())
  if(topo == "L")
    return( line_t_lm())
}

get_lm_time = function(gr,lm = get_topo_lm(gr))
{
  return ( lm$coefficients[1] + size_from_code(gr$code)* lm$coefficients[2])
}

network_from_graph = function( gr)
{

  adj = matrix(rep(0,gr$N*gr$N), ncol = gr$N)

  n_mat = retrace_E(gr ,0)$mat
  link =  ( Mod(n_mat)!= 0 )

  adj[link] = 1

  return ( graph_from_adjacency_matrix(adj))

}

##################################à
# graph search functions

#return distance in link between two sites
distance = function(gr, a = gr$start, b = gr$target)
{
  dist = rep(gr$N+1, gr$N)

  l=1
  dist[a] = 0

  cur_sites = c(a)
  last = 2

  for(i in 1:gr$N)
  {
    #print(cur_sites)

    for(pos in cur_sites)
    {
      conn = ( gr$mat[pos,]!=0 ) & (dist>l )

      #print(dist>l)
      #print(conn)

      dist[conn] = l

      cur_sites = c(cur_sites,(1:gr$N)[conn])
    }

    l= l+1
  }

  #print("Final dist")
  #(dist)

  return( dist[b])
}

sites = function(gr)
{
  return(gr$N)
}

get_list_range = function(gr_list, mode)
{
  if( mode == "sites")
    return( range( sapply(gr_list,sites))  )
  else if( mode == "dist")
    return( range( sapply(gr_list,distance))  )
  else
    print("get_list_range error: mode not found")
}

get_list_x = function(gr_list, mode)
{
  if( mode == "sites")
    return( sapply(gr_list,sites)  )
  else if( mode == "dist")
    return( sapply(gr_list,distance)  )
  else
    print("get_list_x error: mode not found")
}




