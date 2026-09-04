import numpy as np
import matplotlib.pyplot as plt

def exp_expa(x,f,k,alf=1e-3,no_neg=False):
    """
        f=f(x) is the function to be transformed.
        k is a vector of selected values in k-space.
        alf is the strength of regularization. The greater the noise, the higher the alf.
        If we know the coefficients are non-negative, the no_neg is set to True.
    """

    from scipy.optimize import nnls
    from sklearn.linear_model import Ridge

    dk=k[1]-k[0]
    E=np.exp(-x[:,None]*k[None,:])
    E=E*dk
    if no_neg:
        # nk=len(k)
        # reg=np.sqrt(alf)*np.eye(nk)
        # pad=np.zeros(nk)
        # E=np.vstack([E,reg])
        # f=np.concatenate([f,pad])
        coeff,residual=nnls(E,f)
    else:
        ### Ridge Regression / Tikhonov Regularization
        ### Loss=||E*coeff-f||^2+alf*||coeff||^2
        clf=Ridge(alpha=alf,fit_intercept=False)
        clf.fit(E,f)
        coeff=clf.coef_
    frec=np.dot(E,coeff)

    return coeff,frec

dx=0.05
x_len=100
x=np.arange(x_len)*dx
# A=[3.0,1.5]
# K=[2.0,5.0]
A=[3.0,1.5,2.0]
K=[2.0,5.0,8.0]
f=np.sum([A[i]*np.exp(-K[i]*x) for i in range(len(A))],axis=0)
np.random.seed(42)
f=f+0.0*np.random.randn(len(f))
k=np.linspace(0,10,100)
a,f_rec=exp_expa(x,f,k,alf=0.01,no_neg=True)

plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(x,f,'k.',label='Input Data (Noisy)',alpha=0.5)
plt.plot(x,f_rec,'r-',label='Reconstructed f(x)')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(k,a,'b-')
plt.title('Inverted Coefficient Spectrum a(k)')
plt.xlabel('k')
plt.ylabel('Amplitude a(k)')
for i in range(len(K)):
    plt.axvline(K[i],color='g',linestyle='--',label=f'True k={K[i]}')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()




def exp_expa(x,f,k,alf=1e-3,no_neg=False):
    """Small alf: good fit, oscillative solvation, sensitive to noise. Large alf: smooth solvation"""
    from scipy.optimize import lsq_linear
    nk=len(k)

    w=np.zeros(nk)
    w[1:-1]=(k[2:]-k[:-2])/2
    w[0]=(k[1]-k[0])/2
    w[-1]=(k[-1]-k[-2])/2
    E=np.exp(-x[:,None]*k[None,:])
    if alf>0:
        reg=np.eye(nk)*np.sqrt(alf)
        A=np.vstack([E,reg])
        b=np.concatenate([f,np.zeros(nk)])
    else:
        A=E
        b=f
    
    if no_neg:
        res=lsq_linear(A,b,bounds=(0,np.inf),method='trf',tol=1e-10)
    else:
        res=lsq_linear(A,b,bounds=(-np.inf,np.inf),method='trf',tol=1e-10)
    m_sol=res.x
    coeff=m_sol/w
    f_rec=np.dot(E,m_sol)
    
    return coeff,f_rec

dx=0.05
x_len=100
x=np.arange(x_len)*dx
A=[0.3,1.5,10.0]
K=[0.05,0.5,5.0]
# A=[3.0,1.5]
# K=[0.5,5.0]
f=np.sum([A[i]*np.exp(-K[i]*x) for i in range(len(A))],axis=0)
np.random.seed(42)
f=f+0.0*np.random.randn(len(f))
k=np.logspace(-2,1,100)
a,f_rec=exp_expa(x,f,k,alf=1e-5,no_neg=True)

plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(x,f,'k.',label='Input Data (Noisy)',alpha=0.5)
plt.plot(x,f_rec,'r-',label='Reconstructed f(x)')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(k,a,'b-')
plt.title('Inverted Coefficient Spectrum a(k)')
plt.xlabel('k')
plt.ylabel('Amplitude a(k)')
for i in range(len(K)):
    plt.axvline(K[i],color='g',linestyle='--',label=f'True k={K[i]}')
plt.xscale('log')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()