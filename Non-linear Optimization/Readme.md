# Trust Region Method
### Wikipedia definition:
Trust-region method (TRM) is one of the most important numerical optimization methods in solving nonlinear programming (NLP) problems. It works in a way that first define a region around the current best solution, in which a certain model (usually a quadratic model) can to some extent approximate the original objective function. TRM then take a step forward according to the model depicts within the region. Unlike the line search methods, TRM usually determines the step size before the improving direction (or at the same time). If a notable decrease (our following discussion will based on minimization problems) is gained after the step forward, then the model is believed to be a good representation of the original objective function. If the improvement is too subtle or even a negative improvement is gained, then the model is not to be believed as a good representation of the original objective function within that region. The convergence can be ensured that the size of the “trust region” (usually defined by the radius in Euclidean norm) in each iteration would depend on the improvement previously made. 

#### Trust Region Method has been implemented through R; 
following algorithms and steps from the textbook *Optimization – Theory and Practice*. Algorithm of Trust Region Method is outlined as follows: \
<img width="365" alt="image" src="https://user-images.githubusercontent.com/26233980/56931342-add7d380-6aad-11e9-8426-51b0f05b75be.png">

Source code has been scripted according to this algorithm, which has been demonstrated in example five in the textbook. Objective function, used in the midterm, was also tested.
Definition of parameters and functions defined in R:

<img width="471" alt="Table1 1" src="https://user-images.githubusercontent.com/26233980/56931948-da8cea80-6aaf-11e9-8b30-60861b29f7c2.png">
<img width="472" alt="Table1 2" src="https://user-images.githubusercontent.com/26233980/56931956-e11b6200-6aaf-11e9-8933-4e7e392bdcdf.png">

This source code has been specifically coded to follow the same steps as in example five with the objective function given from the midterm. Initial center point and radius was assumed and from this assumptions, interval, gradient, hessian, norm, lower triangular matrix, and empirical threshold values of the ratio  <img width="12" alt="R_k" src="https://user-images.githubusercontent.com/26233980/56932176-ad8d0780-6ab0-11e9-824e-d08eac3c7a93.png">  was calculated for determining the size of the trust-region. 

Norm of the search direction was calculated by “while” loop that is nested inside “if” statement. If the initial norm calculated is not within the demanded interval in this case [0.375, 0.75], then it defines a new function called “new_param”, that calculates new parameter which in turn, calculates new norm until it reaches the demanded interval. 

Threshold value <img width="12" alt="R_k" src="https://user-images.githubusercontent.com/26233980/56932176-ad8d0780-6ab0-11e9-824e-d08eac3c7a93.png"> was calculated with the new search direction and initial center point. Updated lambda, search direction, norm, radius, center point, and demanded interval from the “if-and-while-loop” are displayed under hashtag “new parameters”. 

*Powell’s dogleg method* was also implemented, in this case, the steps it took was much less than it took from original Trust-Region Method. It calculates search direction more directly. As shown in the code, Quasi-Newton direction <img width="16" alt="d_N" src="https://user-images.githubusercontent.com/26233980/56932338-491e7800-6ab1-11e9-8733-51ff457fe0d2.png"> has the same value as the initial search direction. 

Last but not least, minimization of the regularized objective function with its new center point from Trust Region Method, yields the same value as the ‘local model’ from the Trust Region Method with the initial center point.
