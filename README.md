# maxscore-estimator-MIP-inPython

Exact computation of Maximum Score estimator with Mixed Integer Programming via Python code 

Use handy Python code from the file main.py in order to exactly compute the Maxscore estimator with MIP. 
The main script sets up the matrices A,b,c,Aeq,beq,lb,ub of the MIP model, and then relies on a MIP solver to solve it.
We can use CPLEX with the docplex interface. All functions are supplied on top of the single file main.py. 

The dataset is read in readXyw function via the files X.txt, y.txt and w.txt which can be adopted as desired.
Currently, also weights w are supported for an extension called 'maximum weighted score estimator'.
In order to have a more flexible modeling approach, readers are suggested to consult the GAMS version
of the same model in https://www.gams.com/modlib/libhtml/mws.htm.

The code is generic. It supports any reasonable value for sample size, N, and number of predictors, p.

Feedback for the Python code at cflorios@central.ntua.gr, cflorios@aueb.gr.

In case you have trouble using the docplex interface, do not hesitate to contact me for support.

This is a translation of my own Matlab code available in another repository of mine (https://github.com/kflorios/maxscore-estimator-mip).

For completeness, I supply the Matlab manual here too.

The expected result is

value = 765.0
estimates = (1., 170.98512586, 3.65217391, 0.84668192, -110.48970252)
time = 143.0
quality = integer optimal solution


Suggested publication:  

Florios, K., Skouras, S. 
Exact computation of max weighted score estimators
(2008) Journal of Econometrics, 146 (1), pp. 86-91.

http://www.sciencedirect.com/science/article/pii/S0304407608000778 
