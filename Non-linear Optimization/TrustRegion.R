library(pracma)
library(base)
library(MASS)
library(numDeriv)


#Dataset
n <- 10
E <- c(c(9,5),c(6,3),c(9,1),c(8,6),c(8,2),
       c(5,3),c(6,1),c(5,0),c(9,2),c(1,0),
       c(6,5),c(7,6),c(7,5),c(7,1),c(6,0),
       c(4,1),c(3,0),c(7,4),c(9,4),c(9,8),
       c(7,0),c(7,2),c(6,2),c(8,5),c(8,4),
       c(7,3),c(4,0),c(4,2),c(8,0),c(5,1))

#because R, unlike Python, index starts from 1 not 0.
E <- E + 1

y <- matrix(c(0, 0.11, 0, 0.08, 0.06, 0.08, 0.09, 0.04, 0.01, 0,
              0.09, 0, 0, 0, 0.08, 0.05, 0.08, 0.06, 0, 0.03,
              0, 0, 0, 0, 0.11, 0, 0.1, 0.04, 0.04, 0.02,
              0.12, 0, 0, 0, 0, 0.08, 0.09, 0.06, 0, 0,
              0.14, 0.12, 0.09, 0, 0, 0, 0, 0.07, 0.07, 0.03,
              0.12, 0.15, 0, 0.12, 0, 0, 0.07, 0.09, 0.06, 0.01,
              0.11, 0.12, 0.1, 0.11, 0, 0.13, 0, 0.14, 0.05, 0,
              0.16, 0.14, 0.16, 0.14, 0.13, 0.11, 0.06, 0, 0, 0,
              0.19, 0, 0.16, 0, 0.13, 0.14, 0.15, 0, 0, 0.03,
              0, 0.17, 0.18, 0, 0.17, 0.19, 0, 0, 0.17, 0),nrow = 10, ncol = 10, byrow = T)
lambda <- 1

i <- E[seq(1,length(E),2)]
j <- E[seq(2,length(E),2)]
temp <- cbind(i,j)
temp

#objective function
myfunct <- function(theta){
  for (h in 1:nrow(temp)){
    theta_i <- theta[temp[h,1]]              #i-th position
    theta_j <- theta[temp[h,2]]              #j-th position
    L1 <- sum(-y[temp[h,2],temp[h,1]] * (theta_i - theta_j) + log(1+exp(theta_i - theta_j)))
  }
  L2 <- L1 + 0.5 * lambda * norm(theta, "2")^2
  return(L2)
}

#Step 1
#Let center: x^(k) and radius: delta_k, be given;
#declare vars
interval <- list()
center <- list()
radius <- c()
param <- c()
gradi <- list()
Hess <- list()
d <- list()
normi <- c()
L <- list()
w <- list()
w_dotprod <- c()
r_k <- c()


#Trust Region
#Following steps from Example 5 in the textbook.
center[[1]] <- c(1,1,1,1,1,1,1,1,1,1)
radius[1] <- 0.5
param[1] <- 1
  inside <- function(radius){                                                        #setting up the interval
    a <- c(0.75 * radius, 1.5 * radius)
  }
  interval[[1]] <- inside(radius[1])
                                                                                    

minimizer <- function(Hess,param,gradi){                                             #setting up the minimizer
  funct <- -solve(Hess + param * diag(10)) %*% gradi
}


  gradi[[1]] <- grad(myfunct, center[[1]])                                           #grad-update
  Hess[[1]] <- hessian(myfunct, center[[1]])                                         #Hess-update
  
  d[[length(d)+1]] <- minimizer(Hess[[1]], param[length(param)], gradi[[1]])         #minimizer-update
  normi[length(normi)+1] <- norm(d[[length(d)]], "2")                                #norm-update
  
  
  L[[1]] <- t(chol(Hess[[1]]))                                                       #L-update
  w[[1]] <- qr.solve(L[[1]],d[[length(d)]])                                          #w-update
  w_dotprod[1] <- dot(w[[1]],w[[1]])                                                 #dot-update
  
  
#if and while statements to find norm of the search direction(d) to be inside (0.375 , 0.75)
  if ((normi[length(normi)] <= interval[[1]][2] & normi[length(normi)] >= interval[[1]][1]) == FALSE){           
    new_param <- function(param,normi,radius){
      param + (1 - (normi / radius)) * (-(normi)^2 / w_dotprod)
    }
    param[length(param)+1] <- new_param(param[length(param)],normi[length(normi)],radius[1])
    
    d[[length(d)+1]] <- minimizer(Hess[[1]], param[length(param)], gradi[[1]])
    
    normi[length(normi)+1] <- norm(d[[length(d)]], "2")
    
    while ((normi[length(normi)] <= interval[[1]][2] & normi[length(normi)] >= interval[[1]][1]) == F){
      param[length(param)+1] <- new_param(param[length(param)],normi[length(normi)],radius[1])
      
      d[[length(d)+1]] <- minimizer(Hess[[1]],param[length(param)],gradi[[1]])
      
      normi[length(normi)+1] <- norm(d[[length(d)]], "2")
    }
  }
  
  #r_k calc.
  r_k[1] <- (myfunct(center[[1]]) - myfunct(center[[1]]+d[[length(d)]])) / (myfunct(center[[1]]) - (myfunct(center[[1]]) + gradi[[1]] %*% d[[length(d)]] + 0.5 %*% t(d[[length(d)]]) %*% Hess[[1]] %*% d[[length(d)]])) 
  
  if (r_k[1] < 0.25){
    radius[i+1] <- 0.5 * radius[1]
    center[[i+1]] <- center[[1]]
  }else if(r_k[1] > 0.25){
        if(r_k[1] > 0.75 & normi[length(normi)] <= interval[[1]][2] & normi[length(normi)] >= interval[[1]][1]){
          radius[2] <- 2 * radius[1]
        }else{
          radius[i+1] <- radius[1]
        }
    center[[2]] <- center[[1]] + d[[length(d)]]
  }
  interval[[2]] <- inside(radius[2])                                                 #interval-update  

  
#Parameter-Updates displayed.
param
d
normi
radius
center
interval
  
#POWELL's dogleg method.
Pcenter <- list()         #Powells Center Point
#Steepest descent direction d_s
d_s <- -((norm(gradi[[1]], "2")^2) / dot(gradi[[1]],(Hess[[1]] * gradi[[1]]))) * gradi[[1]] 
d_s
x_s <- center[[1]] - d_s    
x_s

#Quasi-Newton direction d_N
d_N <- -solve(Hess[[length(Hess)]]) %*% gradi[[length(gradi)]]
#It has same value as initial search direction from Trust Region Method
d_N
x_N <- center[[1]] + d_N
x_N

#We choose next center point such that:
if (norm(d_N,"2") <= interval[[1]][2] & norm(d_N, "2") >= interval[[1]][1]){
  Pcenter <- x_N
}else if (norm(d_s,"2") >= interval[[1]][2] & norm(d_s, "2") <= interval[[1]][1]){
  Pcenter <- center[[1]] - (radius[1] / norm(gradi[[1]],"2")) %*% gradi[[1]]
}else{
  Pcenter <- (1-param[1])*x_s + param[1]*x_N
}

#new point after dogleg method.
Pcenter

#Testing with Rosenbrock function.
rosenbrock <- function(x){
  n <- length(x)
  x2 <- x[2:n]
  x1 <- x[1:(n-1)]
  sum (100 *(x2 - x1^2)^2 + (x1-1)^2) 
}

rosenbrock(c(1,1,1,1,1,1,1,1,1,1))
rosenbrock(center[[1]])
rosenbrock(center[[2]])

#Minimization of the 'local model' with the initial center point
(myfunct(center[[1]]) + gradi[[1]] %*% d[[length(d)]] + 0.5 %*% t(d[[length(d)]]) %*% Hess[[1]] %*% d[[length(d)]])


#Minimization of the original objective function (midterm function) with the new center point
myfunct(t(center[[2]]))

#They both yield same value







  
  
  
  
  
  