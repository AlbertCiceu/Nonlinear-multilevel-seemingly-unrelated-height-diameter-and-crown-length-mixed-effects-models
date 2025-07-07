fun_HT_CL<-function(dbh,which,a1,a2,b1,b2) {
  value<-which*(1.3 + a1 * (dbh / (1 + dbh)) ^ a2)+
    (1-which)*((dbh / (b1 + exp(b2) * dbh)) ^ 2)
  grad<-cbind(
    a1=which*(dbh/(dbh+1))^a2,
    a2=a1*which*((dbh/(dbh+1))^a2)*log(dbh/(dbh+1)),
    b1=((2*dbh^2)*(which-1))/(exp(b2)*dbh+b1)^3,
    b2=((2*dbh^3)*(which-1)*exp(b2))/(b1+dbh*exp(b2))^3
  )
  attr(value,"gradient") <- grad
  return(value)
}
