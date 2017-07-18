#########################################################################
#               Segmentation Tree: Arbre de segmentation  
#########################################################################

                   # Written by: Nafise Fateri Gouard  #
###########################################################################
# segmentation tree: CART method



# method= la méthode de segmentation
#   .	CHAID : Chi-squared Automatic Interaction Detector
#   .	CART : Classification And Decision Tree
#   .	SIPINA : Système Interactif pour les Processus d'Interrogation Non-Arborescent


Ntree.CART<- function(X, y, xtype, ytype, NaIsModality){
  
  # X : une matrice des variables explicative(here it can be a list of variables, because c(x,y) make every caterorie the same )
  # y: variable à expliqué
  
  # xtype: type des variables explicatives(un vector de)
  # ytype: type de la  variable à expliqué
  
  # Dans la methode CART on utilise différents critère selon le type de Y:
  # le critère de Mesure de l'impureté d'un segment  : 
  #  . Indice de Gini(y nominale)
  #  . Twoingo(y ordinal)
  #  . Least Square Deviation(Moindre carré)(y numeric)
  
  # fixer une règle d'arrêt de la procédure:
  # on arret la procédure quand l'effectif d'un brunch est inférieur à 5% 
  # de l'effectif de la population initial(l'effectif= fréquencie)
  if(NaIsModality==TRUE){for(i in 1:length(X)){X[[i]]<-NaAsModality(X[[i]])}}
  
  if(ytype=="nominal") {fun<-gini} else 
    
  {if(ytype=="ordinal") {fun<-twoing} else {if(ytype=="numeric"){fun<-LSD} 
      
  else {stop("Y's type must be on of: numeric, nominal, ordinal")}}}

  len=0.02*length(na.omit(y))
  continue_tree(X,y,xtype,ytype,fun,parent=0,l=list(),m=0,len)
  
  
}
  
  
    

#--------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------



# i(t) Entropy quadratique
it<-function(y){
  ni=table(y)
  pi=ni/sum(ni)
  P=0
  for(i in 1:length(ni))
  {P=P+pi[i]**2}
  return(1-P)
}



#--------------------------------------------------------------------------------------------------------

# Critère de Gini: y nominal 
gini<-function(x,y){
  #function gini
  # x: un vecteur de la variable(modalité) expliquative
  # y: un vecteur de la variable à expliqué
  
  t=table(x,y)
  S=0
  for(i in 1:dim(t)[1]){
    xx=y[x==names(t[,1])[i]]
    S=S + (sum(t[i,])/sum(t))*it(xx)
  }
  Delta=it(y)-S
  result=list(Delta=Delta)
  return(result)
    
}

#--------------------------------------------------------------------------------------------------------

# critère de LSD, y numeric
LSD<-function(x,y)
{
 # function LSD
 # y<- numeric
 # x<- avec deux modalité ou on le divise juste en deux partie
  
  t=table(x)
  
  tg=y[x==names(t)[1]]
  td=y[x==names(t)[2]]
  ng=sum(table(tg))
  nd=sum(table(td))
  delta=( ng*nd/((ng+nd)**2) ) *( mean(tg,na.rm = TRUE)-mean(td,na.rm = TRUE) )**2
  result=list(Delta=delta)
  return(result)
}
  

#--------------------------------------------------------------------------------------------------------

  
# y ordinal
twoing <- function(x,y){
  # function twoing
  # y: ordinal
  # x: dans deux groupes
  t=table(x,y)
  p=c()
  for( i in 1:length(table(y)) ){
    pg=sum(t[1,1:i])/sum(t[1,])
    pd= sum(t[2,1:i])/sum(t[2,])
    
    pi=(   pg-pd   )**2
    
    p=c(p, pi)
  }
  M=max(p)
  j=which.max(p)
  cut=names(t[1,])[j]
  delta=( sum(t[1,])*sum(t[2,])/(sum(t))**2 )* M 
  result=list(Delta=delta, Cut=cut)
  return(result)
}
  
#------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------

# construire la matrice (dichotomies) de l'ordre m


table.dichotomie<- function(m,s)
 # m= nombre des modalités
 # s= type de variable(numeric, ordinal,nominal)
{ 
  if(!is.integer(m)){m<-floor(m)}
  if(m<2){return(10000)}else{if(m==2){return(nominal=matrix(c(0,1),2, 1))}else{
      if(s=="nominal"){
      t=c(0,1)
      p=2
      
      while(p<m){
        
        z=cbind(t,t,0)
        t=rbind(z,1)
        t[p+1,1:2**(p-1)-1]<-0
        
        p=dim(t)[1]}
        
      }else
        if(s %in% cbind("ordinal", "numeric")){
      ordinal=matrix(0,nrow = m,ncol=m)
      ordinal[lower.tri(ordinal, diag = FALSE)]<-1
      t<-ordinal[,-m]
        }
   else(stop("xtype shoould be one of nominal or ordinal "))
      
      return(t)}
  }
}

#--------------------------------------------------------------------------------------------------------

numeric2ordinal<-function(x,m){
  if(mode(x)!="numeric"){x=as.numeric(x)}
  b<-max(x,na.rm = TRUE)
  a<- min(x,na.rm = TRUE)
  d<-1/(m-2)
  
  for(i in 1:length(x)){ if(!is.na(x[i])) x[i]<-(x[i]-a)/(b-a)}
  
  for(i in 1:length(x)){
    if(!is.na(x[i])) {
      if(x[i]==1){x[i]<-b}else if(x[i]==0){x[i]<-a} else{
      for(k in 1:m-2)
      {if( (k-1)*d < x[i] & x[i]< k*d){x[i]<- a+(k*d-(d/2))*(b-a) }}
    }
    }
  }
  return(x)
}


#---------------------------------------------------------------------------------------------------------

# build brunch(internal calculation)
# this function gives us the best dichotomi with maximum information gain

build_brunches<-function(x,y,xtype,ytype,fun){

  
  if(xtype=="numeric"){x=numeric2ordinal(x,10)}
  t<-table(x)
  m<-length(t) # nombre des modalité de x
  di<-table.dichotomie(m,xtype)
  
  
  n<-names(t)
  delta<-c()
  if(length(dim(di))==0){}else{
  for(k in 1:dim(di)[2] ){
    
    # building 2 groupes
    d<-di[,k]
    tg<-n[d==0]
    td<-n[d==1]
    ind<-matrix(0,1,length(x))
    
    for(jj in 1:length(x))
    { 
      if (is.na(x[jj])){ind[jj]<-NA} else
        if(x[jj] %in% td) {ind[jj]<-1}  
    }
    # fiding la maximume 
    
    f<-fun(ind,y)
    Delta<-f$Delta
    delta=c(delta, Delta)  
  }
  
  informationgain=max(delta)
  best.dichotomy.index<-which.max(delta)
  d<-di[,best.dichotomy.index]
  tg<-n[d==0]
  td<-n[d==1]
  { ind<-matrix(0,1,length(x))
    for(jj in 1:length(x))
      if (is.na(x[jj])){ind[jj]<-NA} else
        if(x[jj] %in% td) {ind[jj]<-1}  
  }

  result=list(new_x=ind,informationgain=informationgain,tg=tg,td=td)
  return(result)
  }
}  

#-------------------------------------------------------------------------------------------------------
best_brunch<-function(X,y,xtype,ytype,fun){
# X est matrice des variables explicatives
# xtype unecolonnes des types des variables explicatives
info=c()
for(j in 1:length(X)){
  # x est une variable explicative
  # xtype<- type de x : nominal, ordinal
  brunches<-build_brunches(X[[j]],y,xtype[j],ytype,fun)
  #compare the information gains and take the maximum one
  info=cbind(info, brunches$informationgain )
}
Delta=max(info)
best.dichotomy.variable<-which.max(info)

best_brunch<-build_brunches(X[[best.dichotomy.variable]],y,xtype[best.dichotomy.variable],ytype,fun)
result<-list(all=info,best_variable=names(X[best.dichotomy.variable]),brunch=best_brunch)
return(result)
}

$brunch$informationgain< 0.0001

#------------------------cycle functoin-----------------------



continue_tree<-function(X,y,xtype,ytype,fun,parent,l=list(),m=0,len){
  print(m)
  print('======================================================================================================')
  
  m<-m+1
  print(m)
  print(inter_calcul(y, ytype))
  print(table(y))
  print('-----------------------------------------------------')
 
  if(parent<length(X)+1 & length(na.omit(y)!=0) & length(y)> len  & all(table(y)>2)){
      parent <-parent+1
      nn <-inter_calcul(y, ytype)
      #print(nn)
      l[[parent]]<-list()
      l[[parent]][[1]]<-nn
        new_brunch <- best_brunch(X,y,xtype,ytype,fun)
       if( !is.null (new_brunch$brunch$informationgain) ){
          resultatg=list(new_brunch$best_variable,Delta=new_brunch$brunch$informationgain ,g=new_brunch$brunch$tg)
          print(resultatg)
          resultatd=list(new_brunch$best_variable,Delta=new_brunch$brunch$informationgain ,d=new_brunch$brunch$td)
          print(resultatd)
          print('************************')
           #l[[parent]][[1]]<-resultatg
           #l[[parent]][[2]]<-resultatd
           # (l)

           #print(resultatg)
           #print(resultatd)
        
        
          #if(length(nn$ni)==0){print("gouard jan finish")}else{
          #  if(sum(nn$ni)<5 | length(unique(na.omit(y)))==1){print("nafis joonam khoobe")} else
          # in this step, we need to remove the best variable's dichotomy from the matrix of variables
 
           xx<-new_brunch$brunch$new_x
           new_y<-split(y,xx)
           Xg<-list()
           Xd<-list()
           j=0;
           for(i in 1:length(X)){
               j=j+1;
               new_x<-split(X[[i]],xx)
               Xg[[i]]<-new_x$'0'
               Xd[[i]]<-new_x$'1'
               names(Xg)[i]<-names(X)[i]
               names(Xd)[i]<-names(X)[i]
           }
           
  
      #a<-paste(a,'g',parent)
      continue_tree(Xg,new_y$`0`,xtype,ytype,fun,parent,l,m,len)
      #b<-paste(b,'d',parent)
      continue_tree(Xd,new_y$`1`,xtype,ytype,fun,parent,l,m,len)
 
  
  } }
return(l)
}


#---------------inter_calcule--------------------------------------

inter_calcul<-function(y,ytype){
  if(ytype=='nominal') {return(y=table(y))} 
  else{
  m<-mean(y,na.rm=TRUE)
  y<-yOrdered(y)
  t<-table(y)
  ni=t
  pi=round(ni/sum(t),digits=2)
  NPS<-pi[names(pi)=="p"]-pi[names(pi)=="d"]
  return(list(NPS=as.numeric(NPS),CSAT=m))
}

}


















############################## Testing the data #########################
yOrdered<-function(x){
  for(i in 1:length(x)){
    if(is.na(x[i])){}
    else if(x[i]==9 | x[i]==10) {x[i] <-"p"} else if(x[i]<7) {x[i] <-"d"} else x[i]<-"n"
  }
  return(x)
  
}

#********************exemple: Air France******************************

 xtype=cbind("nominal","ordinal","ordinal","ordinal")
 ytype="ordinal"
 y=data3$Q5_NPS
 #y<-yOrdered(y)
 x1<-data3$Q4_FCR
 x2<-data3$Q3_CLARTE
 x3<-data3$Q3_AMABILITE
 #x4<-data3$Q1_SATISFACTION_GLOBALE
 X=list()
 X[[1]]<-x1
 X[[2]]<-x2
 X[[3]]<-x3
 X[[4]]<-x4
 name<-c("FCR","Clarté","amabilité","satis_globale")
 names(X)<-name
 
 bb<-Ntree.CART(X,y,xtype,ytype,NaIsModality=TRUE)
 
 
 
 library(rpart)
 library(rpart.plot)
 library(ggplot2)
 m2<- rpart(Q5_NPS~Q3_CLARTE+Q4_FCR+Q3_AMABILITE, data = data3 ,method="class")
 #m2<- rpart(Q5_NPS~Q3_CLARTE+Q4_FCR+Q3_AMABILITE+ Q1_SATISFACTION_GLOBALE, data = data3 ,method="class")
 rpart.plot(m2, type=3, digits =  3 , fallen.leaves = TRUE )
 

 
#**************************exemple: idTGV ************************** 
 
xtype=cbind("numeric","nominal","nominal")
ytype="ordinal"
y=data2$RECOMMANDATION
x1=data2$retards
x2=data2$site_reservation
x3=data2$ambiance
X=list()
X[[1]]<-x1
X[[2]]<-x2
X[[3]]<-x3



name<-c("retard","site_reservation","ambiance")
names(X)<-name


bb<-Ntree.CART(X,y,xtype,ytype,NaIsModality=TRUE)


m2<- rpart(RECOMMANDATION~retards+site_reservation+ambiance  , data =idTGV  ,method="class")
rpart.plot(m2, type=3, digits =  3 , fallen.leaves = TRUE )



####################################
x
#"""""""""""""""""""""""""""""""""""""""""""
y<- kyphosis$Kyphosis
x1<-kyphosis$Age
x2<-kyphosis$Number
x3<-kyphosis$Start


X=list()
xtype=cbind("numeric","numeric","numeric")
X[[1]]<-x1
X[[2]]<-x2
X[[3]]<-x3
ytype="nominal"
names(X)<-cbind("Age","number","start")

bb<-Ntree.CART(X,y,xtype,ytype,NaIsModality = TRUE)



m2<- rpart(Kyphosis~Age+Number+Start  , data =kyphosis  ,method="class")
rpart.plot(m2, type=3, digits =  3 , fallen.leaves = TRUE )
plot(y~x3)
######################################################


NaAsModality<-function(X)
{
 X[is.na(X)]<-"NA"
 return(X)
}






###########################################################


tet=list()

test.tet <- function() {
  6+5
}
for(i in 1:5){
  tet[[i]]<-i+10
  
  #new_x<-split(X[[i]],xx)
  #Xg[[i]]<-new_x$'0'
  #Xd[[i]]<-new_x$'1'
  #names(Xg)[i]<-names(X)[i]
  #names(Xd)[i] <-names(X)[i]
}
