# maxscore-estimator-MIP-inPython

Exact computation of Maximum Score estimator with Mixed Integer Programming via Python code 

Use handy Python code from the file main.py in order to exactly compute the Maxscore estimator with MIP. 
The main script sets up the matrices A,b,c,Aeq,beq,lb,ub of the MIP model, and then relies on a MIP solver to solve it.
We can use CPLEX with the docplex interface. All functions are supplied on top of the single file main.py. 

The dataset is read in readXyw function via the files X.txt, y.txt and w.txt which can be adopted as desired.
Currently, also weights w are supported for an extension called 'maximum weighted score estimator'.
In order to have a more flexible modeling approach, readers are suggested to consult the GAMS version
of the same model in https://www.gams.com/modlib/libhtml/mws.htm.

The code is semi-generic. It supports arbitrary values for sample size, N, and number of predictors, p.

The only piece of code the user has to edit is these two lines of code in milp_cplex function:

    for i in R1:
        model.add_constraint((XX[i,0]*beta[0]+XX[i,1]*beta[1]+
                                  XX[i,2]*beta[2]+XX[i,3]*beta[3]+
                                  XX[i,4]*beta[4]) + M[i]*z[i]  <= M[i] )

Where I have hard coded the case with p=5.

Feel free to edit these two lines, and facilitate any number of predictors p. Usually p is a small integer (2 to say 10).

Feedback for the Python code at cflorios@central.ntua.gr, cflorios@aueb.gr.

In case you have trouble using the docplex interface, do not hesitate to contact me for support.

This is a translation of my own Matlab code available in another repository of mine (https://github.com/kflorios/maxscore-estimator-mip).

For completeness, I supply the Matlab manual here too.

The expected result is

value = 764.9999999999998

estimates = (1., 170.98512586, 3.65217391, 0.84668192, -110.48970252)

time = 128.43699999999808

quality = integer optimal solution


Suggested publication:  

Florios, K., Skouras, S. 
Exact computation of max weighted score estimators
(2008) Journal of Econometrics, 146 (1), pp. 86-91.

http://www.sciencedirect.com/science/article/pii/S0304407608000778 
