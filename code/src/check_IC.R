check_IC<-function(mat_inf, mat_sup, mat_true, mat_bin){
S<-dim(mat_bin)[1]
r<-dim(mat_bin)[2]

  for(i in 1:S){
    for(j in 1:r){
      if(mat_inf[i,j]<=mat_true[i,j] && mat_sup[i,j]>=mat_true[i,j]){
        mat_bin[i,j]=1
      }
    }
  }
  mat_bin
}


