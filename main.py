import cplex
import sys
import numpy as np
import pandas as pd
import math
def MaxScoreCompute():
    #Main Program: Computes Max score by defining a MILP and calling milp.py
    X,y,w = readXyw()
    X,mu,sigma = standardizeX(X)
    c,A,b = definecAb(X,y,w)
    lb,ub,Aeq,beq,n,p,best = definelbub(X)
    x,score,feasible,time = milp_cplex(c,A,b,Aeq,beq,lb,ub,n,p,X,y, best)
    estimatesNorm=x
    value=score
    status=feasible
    runtime=time
    quality=status
    estimatesRaw=denormalizeEstimates(estimatesNorm,mu,sigma)
    estimates=estimatesRaw
    return value,estimates,time,quality


def readXyw():
    #Reads X,y,w of given max score problem
    X=np.loadtxt("X_Horowitz.txt")
    #X = np.loadtxt("X_numeric200.txt")
    X=X[:,1:]
    y=np.loadtxt("y_Horowitz.txt")
    #y = np.loadtxt("y_numeric200.txt")
    y=y[:,1:]
    w=np.loadtxt("w_Horowitz.txt")
    #w = np.loadtxt("w_numeric200.txt")
    w=w[:,1:]
    return X,y,w

def standardizeX(X):
    #Standardizes X
    df = pd.DataFrame.from_records(X)
    mu = df.mean(axis=0)
    sigma = df.std(axis=0)
    testX = (X - np.tile(mu,(np.size(X,0),1))) / np.tile(sigma,(np.size(X,0),1))
    p = np.size(X,1)
    for j in range(p):
        if math.isnan(testX[0,j]):
            X[:,j] = X[:,j]
        else:
            X[:,j] = testX[:,j]
    return X,mu,sigma

def definecAb(X,y,w):
    #Defines c,A,b for milp.py
    n=np.size(X,0)
    p=np.size(X,1)
    c1=np.tile(-1,(1,n))
    c2=np.tile(0,(1,p))
    c=np.hstack((c1,c2))
    #d=15
    d=10
    #d=5
    M=np.zeros(n)
    for i in range(n):
        #M[i]=np.abs(X[i,0])+np.abs(X[i,1:]@np.transpose(np.tile(d,(1,p-1))))
        #M[i] = np.abs(X[i, 0]) + np.abs(X[i, 1])* d + np.abs(X[i, 2])* d +np.abs(X[i, 3])* d +np.abs(X[i, 4])* d
        #M[i] = np.abs(X[i, 0]) + np.sum(np.abs(X[i,j])*d for j in range(1,p))
        M[i] = np.abs(X[i, 0]) + sum(np.abs(X[i, j]) * d for j in range(1, p))
    Abin=np.diag(M)
    Areal=np.zeros((n,p))
    for i in range(n):
        for j in range(p):
            Areal[i,j]=(1-2*y[i])*X[i,j]
    A = np.concatenate((Abin,Areal),axis=1)
    b = M
    return c,A,b

def definelbub(X):
    #Defines lb,ub for milp.py
    #d=15
    d=10
    #d=5
    n=np.size(X,0)
    p=np.size(X,1)
    lb1=np.tile(0,(1,n))
    lb2=np.tile(-d,(1,p))
    lb2[0,0]=1
    lb=np.hstack((lb1,lb2))
    ub1=np.tile(1,(1,n))
    ub2=np.tile(d,(1,p))
    ub2[0,0]=1
    ub=np.hstack((ub1,ub2))
    Aeq=None
    beq=None
    best=0
    return lb,ub,Aeq,beq,n,p,best

def milp_cplex(c,A,b,Aeq,beq,lb,ub,n,p, X,y,best):
    # Solves a mixed integer lp using cplex 20.1
    # c: is objective function coefficients A: is constraint matrix
    # b: is constraint vector
    # lb: lower bound ub: upper bound n: number of 0-1 variables
    # best: is best solution so far
    # Note this uses the Python/Cplex interface documented at
    # https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html
    # Also, it assumes first n variables must be integer.
    # The MIP equations of the maximum score estimator are available at
    # Florios. K, Skouras, S. (2008) Exact computation of maximum weighted
    # score estimators, Journal of Econometrics 146, 86-91.
    # Written by Kostas Florios, December 26, 2020
    #
    #
    # echo %PYTHONPATH%
    # C:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\python\3.7\x64_win64
    # (setup a conda environment called "cplex")
    # conda create --name=cplex python=3.7.6
    # conda activate cplex
    # run the proper install procedure in that environment (see ibm link above)
    # also pip install numpy and pip install pandas in the conda environment
    # python main.py

    from docplex.mp.model import Model
    M = b

    model=Model("maxscore")
    R1 = range(n)
    R2 = range(p)
    idx1  = [(i) for i in R1]
    idx2 = [(j) for j in R2]

    #d=15
    d=10
    #d=5
    z = model.binary_var_dict(idx1, None)
    beta = model.continuous_var_dict(idx2,lb=-d,ub=d)

    XX =np.zeros((n,p))
    for i in R1:
        for j in R2:
            XX[i,j]=(1-2*y[i])*X[i,j]

    for i in R1:
        #model.add_constraint((XX[i,0]*beta[0]+XX[i,1]*beta[1]+
        #                          XX[i,2]*beta[2]+XX[i,3]*beta[3]+
        #                          XX[i,4]*beta[4]) + M[i]*z[i]  <= M[i] )
        model.add_constraint(M[i]*z[i]+model.sum(XX[i,j]*beta[j] for j in R2) <= M[i] )


    #this doesn't work
    #for i in R1:
    #    if (y[i] ==0):
    #        model.add_constraint((X[i,0]*beta[0]+X[i,1]*beta[1]+
    #                              X[i,2]*beta[2]+X[i,3]*beta[3]+
    #                              X[i,4]*beta[4]) + M[i]*z[i]  <= M[i] )
    #    if (y[i] == 1):
    #        model.add_constraint((X[i, 0] * beta[0] + X[i, 1] * beta[1] +
    #                              X[i, 2] * beta[2] + X[i, 3] * beta[3] +
    #                              X[i, 4] * beta[4]) + M[i] * z[i] >= M[i])

    model.add_constraint(beta[0]==1)
    for j in R2:
        model.add_constraint(beta[j] >= -d)
        model.add_constraint(beta[j] <= +d)

    model.total_inside_obj = model.sum(z[i] for i in R1)
    model.add_kpi(model.total_inside_obj, "inside cost")
    model.maximize(model.total_inside_obj)
    model.print_information()
    model.export_as_lp("model.lp")
    ok=model.solve()
    print(ok)
    #x=ok.get_value_list([beta[0],beta[1],beta[2],beta[3],beta[4]])
    x = ok.get_value_list(list(beta[j] for j in range(p)))
    score=ok.objective_value
    feasible=ok.solve_details.status
    time=ok.solve_details.time
    return x,score,feasible,time

def denormalizeEstimates(estimatesNorm,mu,sigma):
    #denormalized estimatesNorm obtained by Cplex MIP to estimatesRaw
    #which are meaningful to the user

    #quick and dirty implementation, based on GAMS and Fortran Analogues
    p = len(estimatesNorm)
    betaNorm=estimatesNorm
    betaRaw = np.zeros(p)
    betaHelp = np.zeros(p)
    for j in range(p):
        if (sigma[j] != 0):
            betaHelp[j] = betaNorm[j]/sigma[j]
        if (sigma[j] ==0):
            for jj in range(p):
                if (sigma[jj] != 0):
                    betaHelp[j] = betaHelp[j] - betaNorm[jj]*mu[jj]/sigma[jj]
                else:
                    jj0=jj
            betaHelp[j] = betaHelp[j]+betaNorm[jj0]

    for j in range(p):
        betaRaw[j] = betaHelp[j]/betaHelp[0]

    estimatesRaw=betaRaw
    return estimatesRaw


if __name__ == "__main__":

    #MaxScoreCompute()
    value, estimates, time, quality = MaxScoreCompute()
    print(value)
    print(estimates)
    print(time)
    print(quality)