# Arbre de regression à ma facon



function(x,y,n){
  # x<---c'est une matrice des variables
  # y<---c'est le key-variable non codé
  # n<-- dans quelle étap il faut arreter?
  
  yy<-ord2detpro(y)
  x$NPS<-yy
  t<-table.cor(x,y)
  z=1
  
  p=dim(x)[2]
  
  zStore=list()
    
  for(i in 1:p){
     
    j=as.numeric(t[i+1,3])
    
    
    for(k in 1:length(z)){
      
      if(i==1){ # premier étape
        
        Inter.calcul(x,x[,p],j)
        zStore[[i]]=split(x,x[,j])
        # print(names(zStore[[1]]))
      }
      
      else{ # les étapes suivantes
        
        
        xx<-as.data.frame(zStore[k])
        
      
      
      Inter.calcul(x,x[,p],j)

      zStore[k]<- split(x,x[,j]) # c'est une class inclut les data.frames (length(z)=modalite(x,j))
      }
      
    }
    
  }
}

classify=function(x,j){ #x est une form d'un list
  
  for(k in 1:length(x))
    x=as.matrix(x[k])
    split(x, x[,j])
    Inter.calcul(x,j)
    
    print
}

Inter.calcul<-function(x,j){#le NPS est toujours à la fin
  p=dim(x)[2]
  tt=table(x[,c(j,p)]) # inja moshkel vojoood dare nemidoonam chera bayad vaghty se shanbe oomadi in barname ro dorost koni
  a=apply(tt,2,sum)/sum(tt) # le pourcentage de Prpmoteur et detrecteur en root
  b=apply(tt,1,sum)/sum(tt) # le pourcentage sur les branches
  
  print(a)
  print(b)
  
}

