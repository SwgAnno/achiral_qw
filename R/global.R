# a bunch of those methods depend on the global variable N which
# stores the number of sites of the ring only

# N = 7
# phi = 1


#this scripd is no longer to be referenced apar from the setup function
# we no longer use a global variable H

#ring = function(N)
#{
#  mat = matrix(complex(N*N), ncol = N)
#
#  for( i in 1:N)
#    {
#    mat[i,i] = as.complex(2)
#    mat[i,(i%%N+1)] <- as.complex(-1)
#    mat[(i%%N+1),i] <- as.complex(-1)
#    }
#
#  assign("N", N, envir = .GlobalEnv)
#  mat
#}

#ringNhandle = function(N)
#{
#  mat = matrix(complex((N+2)*(N+2)), ncol = N+2)
#
#  for( i in 1:N)
#    {
#    mat[i,i] = as.complex(2)
#    mat[i,(i%%N+1)] <- as.complex(-1)
#    mat[(i%%N+1),i] <- as.complex(-1)
#    }
#
#  mat[N+1,N+1] = 2
#  mat[N+2,N+2] = 2
#
#  mat[N+1,1] = -1
#  mat[1,N+1] = -1
#
#  mat[N+2,N%/%2+1] = -1
#  mat[N%/%2+1,N+2] = -1
#
#  assign("N", N, envir = .GlobalEnv)
#
#  mat
#}
#
#rephase = function( phi = exp(1.57i))
#{
#
#  m[1,N] = -1*phi
#  m[N,1] = -1*Conj(phi)
#
#  m
#}

nsub = function(N,codes)
  return(  gsub("N",N, codes))


prepare_globals = function()
{

	assign("H", m, envir = .GlobalEnv)
	assign("H_ev", eigen(H), envir = .GlobalEnv)

}

setup = function()
{
	library(rootSolve)
	library(stringr)
	library(igraph)
  library(TeachingDemos)
  library(fields)
  library(assertive.types)

	src = c( "global.R", "evoluzione.R", "objMat.R", "minmax.R", "minmax_2.R", "simulaMat.R", "lm.R", "varie.R", "line_sub.R")
	assign("files", src, envir = .GlobalEnv)

	assign("GLOBAL_TIME_SEARCH_CONSTANT", 1, envir = .GlobalEnv)
  assign("GLOBAL_FIRST_MAX", TRUE, envir = .GlobalEnv)
	assign("GLOBAL_UNIROOT_TOL",.Machine$double.eps^0.2, envir = .GlobalEnv)
  assign("GLOBAL_RANDOM_START_COUNT",30, envir = .GlobalEnv)

  assign("YOLO", FALSE, envir = .GlobalEnv)

  assign("GLOBAL_RANDOM_START_COUNT",30, envir = .GlobalEnv)

  assign("p_max_func", p_max_first, envir = .GlobalEnv)
  assign("t_p_max_func", t_p_max_first, envir = .GlobalEnv)

  assign("globalUtils",new.env(), envir = .GlobalEnv)

  #TC stack
  globalUtils$TCstack = rep(1,10)
  globalUtils$TCpoint = 1


  #FC stack
  globalUtils$FCstack = rep(TRUE,10)
  globalUtils$FCpoint = 1

  #FIRST stack

	print(paste("Tolleranza algoritmo di ricerca = ", GLOBAL_UNIROOT_TOL))

}

update_source = function()
{
	sapply(files,source)
}

################################
## TCstack utilities

getTC = function()
{
  return(globalUtils$TCstack[globalUtils$TCpoint])
}

setTC = function(val)
{
  globalUtils$TCstack[globalUtils$TCpoint] = val
}

pushTC = function(val)
{
  globalUtils$TCpoint = globalUtils$TCpoint+1
  globalUtils$TCstack[globalUtils$TCpoint] = val
}

popTC = function()
{
  if( globalUtils$TCpoint != 1)
  {
    globalUtils$TCpoint = globalUtils$TCpoint-1
  }

  else
  {
    print("TCstack ERROR: Trying to pop 1 item stack")
  }
}

################################
## FCstack utilities

getFC = function()
{
  return(globalUtils$FCstack[globalUtils$FCpoint])
}

setFC = function(val)
{
  globalUtils$FCstack[globalUtils$FCpoint] = val

  if(getFC())
  {
    assign("p_max_func", p_max_first, envir = .GlobalEnv)
    assign("t_p_max_func", t_p_max_first, envir = .GlobalEnv)
  }
  else
  {
    assign("p_max_func", p_max, envir = .GlobalEnv)
    assign("t_p_max_func", t_p_max, envir = .GlobalEnv)
  }
}

switchFC = function()
{
  setFC(!globalUtils$FCstack[globalUtils$FCpoint])
}

pushFC = function(val)
{
  globalUtils$FCpoint = globalUtils$FCpoint+1
  setFC(val)
}

popFC = function()
{
  if( globalUtils$FCpoint != 1)
  {
    globalUtils$FCpoint = globalUtils$FCpoint-1

    setFC(getFC())
  }

  else
  {
    print("FCstack ERROR: Trying to pop 1 item stack")
  }
}

################################

