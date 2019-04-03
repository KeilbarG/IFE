
tensor <- function(AA)
{
nr = nrow(AA)
nl = ncol(AA)
phi_tensor = matrix(0, nr^2, nl^2)
for(i in 1:nr){
  for(j in 1:nl){
    phi_tensor[(nl*(i-1)+1):(nr*i),(nl*(j-1)+1):(nr*j)] = AA[i,j] * AA
  }
} 
return(phi_tensor)

}

vec  <- function(AA)

{
 

 nr = nrow(AA)
 nl = ncol(AA)
 eps_vec = matrix(0, nr*nl, 1)
 for(i in 1:nr){
     eps_vec[(nr*(i-1)+1):(nr*i)] = Sigma_eps[,i]
    }
  return(eps_vec)

}