                         ## controlled
                         ########SATELLITE#######
                            ##T=0.1
import time
import numpy as np
import shamirproc
import shamirscheme as sss
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from tabulate import tabulate
from decimal import Decimal
from scipy.io import savemat

def control(r,Tf,F,server,t,scale,sigma,stepinfo):
    


    #initial values
    x1=0
    x2,x3,x4,y=0,0,0,0
    X1=[]
    X2=[]
    X3=[]
    X4=[]
    Y=[]
    Tc=1
    dt=0.005
    samp=Tf/dt
    samp=int(samp)
    samples=np.linspace(0,Tf, num=samp)

    #noise
    noise=np.random.normal(0,sigma,samp)

    #UNCONTROLLED SYSTEM
    for i in range(samp):
        x1=0.9998*x1-0.004999*x2+0.004999*Tc
        x2=0.004999*x1+1*x2+1.25e-05*Tc
        x3=1.25e-05*x1+0.005*x2+1*x3+2.083e-08*Tc
        x4=2.083e-08*x1+1.25e-05*x2+0.005*x3+1*x4+2.604e-11*Tc
        y=0.036*x3+0.9*x4
        X1.append(x1)
        X2.append(x2)
        X3.append(x3)
        X4.append(x4)
        Y.append(y)




    #controlled
    #initial values
    x1c,x2c,x3c,x4c,yc,u=0,0,0,0,0,0
    X1c=[]
    X2c=[]
    X3c=[]
    X4c=[]
    Yc=[]
    U=[]
    z1c,z2c,z3c,z4c,uc=0,0,0,0,0
    Z1c=[]
    Z2c=[]
    Z3c=[]
    Z4c=[]
    Uc=[]
    kr=0.3933
    tic=time.time()
    for i in range(samp):

        u=r*kr-uc
        x1c=0.9998*x1c-0.004999*x2c+0.004999*u
        x2c=0.004999*x1c+1*x2c+1.25e-05*u
        x3c=1.25e-05*x1c+0.005*x2c+1*x3c+2.083e-08*u
        x4c=2.083e-08*x1c+1.25e-05*x2c+0.005*x3c+1*x4c+2.604e-11*u
        yc=0.036*x3c+0.9*x4c

        z1c=1*z1c-1.27e-05*z2c-0.0009439*z3c-0.02316*z4c+0.0001303*u-0.006641*yc
        z2c=0.04*z1c+1*z2c-0.002388*z3c-0.05874*z4c+0.0001813*u-0.01312*yc
        z3c=0.0007996*z1c+0.03996*z2c+0.997*z3c-0.0733*z4c+0.0001282*u-0.003186*yc
        z4c=2.075e-05*z1c+0.001542*z2c+0.07567*z3c+0.8928*z4c+9.078e-05*u+0.5447*yc
        uc=64*z4c

        X1c.append(x1c)
        X2c.append(x2c)
        X3c.append(x3c)
        X4c.append(x4c)
        Yc.append(yc)
        Z1c.append(z1c)
        Z2c.append(z2c)
        Z3c.append(z3c)
        Z4c.append(z4c)
        U.append(u)
        Uc.append(uc)
    toc=time.time()

    ###########################################noisy controlled#######################################################
    x1cn,x2cn,x3cn,x4cn,ycn,un=0,0,0,0,0,0
    X1cn=[]
    X2cn=[]
    X3cn=[]
    X4cn=[]
    Ycn=[]
    Un=[]
    z1cn,z2cn,z3cn,z4cn,ucn=0,0,0,0,0
    Z1cn=[]
    Z2cn=[]
    Z3cn=[]
    Z4cn=[]
    Ucn=[]
    kr=0.3933

    for i in range(samp):

        un=r*kr-ucn
        x1cn=0.9998*x1cn-0.004999*x2cn+0.004999*un
        x2cn=0.004999*x1cn+1*x2cn+1.25e-05*un
        x3cn=1.25e-05*x1cn+0.005*x2cn+1*x3cn+2.083e-08*un
        x4cn=2.083e-08*x1cn+1.25e-05*x2cn+0.005*x3cn+1*x4cn+2.604e-11*un
        ycn=0.036*x3cn+0.9*x4cn+noise[i]

        z1cn=1*z1cn-1.27e-05*z2cn-0.0009439*z3cn-0.02316*z4cn+0.0001303*un-0.006641*ycn
        z2cn=0.04*z1cn+1*z2cn-0.002388*z3cn-0.05874*z4cn+0.0001813*un-0.01312*ycn
        z3cn=0.0007996*z1cn+0.03996*z2cn+0.997*z3cn-0.0733*z4cn+0.0001282*un-0.003186*ycn
        z4cn=2.075e-05*z1cn+0.001542*z2cn+0.07567*z3cn+0.8928*z4cn+9.078e-05*un+0.5447*ycn
        ucn=64*z4cn

        X1cn.append(x1cn)
        X2cn.append(x2cn)
        X3cn.append(x3cn)
        X4cn.append(x4cn)
        Ycn.append(ycn)
        Z1cn.append(z1cn)
        Z2cn.append(z2cn)
        Z3cn.append(z3cn)
        Z4cn.append(z4cn)
        Un.append(un)
        Ucn.append(ucn)    
    
    #############################################secure control############################################################
    rs=sss.share(F,r,t,server,scale)
    x1cs,x2cs,x3cs,x4cs,ycs,us,us_rec=0,0,0,0,0,0,0
    X1cs=[]
    X2cs=[]
    X3cs=[]
    X4cs=[]
    Ycs=[]
    Us=[]
    Us_rec=[]

    ucs=sss.share(F,0,t,server,scale)
    z1cs=sss.share(F,0,t,server,scale)
    z2cs=sss.share(F,0,t,server,scale)
    z3cs=sss.share(F,0,t,server,scale)
    z4cs=sss.share(F,0,t,server,scale)

    Z1cs=[]
    Z2cs=[]
    Z3cs=[]
    Z4cs=[]
    Ucs=[]
    Z1cs_rec=[]
    Z2cs_rec=[]
    Z3cs_rec=[]
    Z4cs_rec=[]
    Ucs_rec=[]


    #sss.R2F(F,5.2323,scale)

    Ac11=sss.R2F(F,1,scale)
    Ac12=sss.R2F(F,-1.27e-05,scale)
    Ac13=sss.R2F(F,-0.0009439,scale)
    Ac14=sss.R2F(F,-0.02316,scale)

    Ac21=sss.R2F(F,0.04,scale)
    Ac22=sss.R2F(F,1,scale)
    Ac23=sss.R2F(F,-0.002388,scale)
    Ac24=sss.R2F(F,-0.05874,scale)

    Ac31=sss.R2F(F,0.0007996,scale)
    Ac32=sss.R2F(F,+0.03996,scale)
    Ac33=sss.R2F(F,0.997,scale)
    Ac34=sss.R2F(F,-0.0733,scale)

    Ac41=sss.R2F(F,2.075e-05,scale)
    Ac42=sss.R2F(F,0.001542,scale)
    Ac43=sss.R2F(F,0.07567,scale)
    Ac44=sss.R2F(F,0.8928,scale)

    Bc11=sss.R2F(F,+0.0001303,scale)    
    Bc12=sss.R2F(F,-0.006641,scale)

    Bc21=sss.R2F(F,+0.0001813,scale)
    Bc22=sss.R2F(F,-0.01312,scale)

    Bc31=sss.R2F(F,+0.0001282,scale)
    Bc32=sss.R2F(F,-0.003186,scale)

    Bc41=sss.R2F(F,+9.078e-05,scale)
    Bc42=sss.R2F(F,+0.5447,scale)

    Cc14=sss.R2F(F,64,scale)

    #Ac11=sss.share(F,1,t,server,scale)
    #Ac12=sss.share(F,-1.27e-05,t,server,scale)
    #Ac13=sss.share(F,-0.0009439,t,server,scale)
    #Ac14=sss.share(F,-0.02316,t,server,scale)

    #Ac21=sss.share(F,0.04,t,server,scale)
    #Ac22=sss.share(F,1,t,server,scale)
    #Ac23=sss.share(F,-0.002388,t,server,scale)
    #Ac24=sss.share(F,-0.05874,t,server,scale)

    #Ac31=sss.share(F,0.0007996,t,server,scale)
    #Ac32=sss.share(F,+0.03996,t,server,scale)
    #Ac33=sss.share(F,0.997,t,server,scale)
    #Ac34=sss.share(F,-0.0733,t,server,scale)

    #Ac41=sss.share(F,2.075e-05,t,server,scale)
    #Ac42=sss.share(F,0.001542,t,server,scale)
    #Ac43=sss.share(F,0.07567,t,server,scale)
    #Ac44=sss.share(F,0.8928,t,server,scale)

    #Bc11=sss.share(F,-0.0001303,t,server,scale)   
    #Bc12=sss.share(F,+0.006641,t,server,scale)

    #Bc21=sss.share(F,-0.0001813,t,server,scale)
    #Bc22=sss.share(F,+0.01312,t,server,scale)

    #Bc31=sss.share(F,-0.0001282,t,server,scale)
    #Bc32=sss.share(F,+0.003186,t,server,scale)

    #Bc41=sss.share(F,-9.078e-05,t,server,scale)
    #Bc42=sss.share(F,-0.5447,t,server,scale)

    #Cc14=sss.share(F,64,t,server,scale)



    krs=sss.share(F,0.3933,t,server,scale)
    rec=sss.basispoly(F,server)
    
    tics=time.time()
    for i in range(samp):
        us=shamirproc.add(F,shamirproc.mult(F,rs,krs,t,server,scale),-ucs,t,server,scale)
        Us.append(us)
        us_rec=sss.rec(F,us,rec,scale)
        Us_rec.append(us_rec)
        x1cs=0.9998*x1cs-0.004999*x2cs+0.004999*us_rec
        x2cs=0.004999*x1cs+1*x2cs+1.25e-05*us_rec
        x3cs=1.25e-05*x1cs+0.005*x2cs+1*x3cs+2.083e-08*us_rec
        x4cs=2.083e-08*x1cs+1.25e-05*x2cs+0.005*x3cs+1*x4cs+2.604e-11*us_rec
        ycs=0.036*x3cs+0.9*x4cs
        Ycs.append(ycs)
        X1cs.append(x1cs)
        X2cs.append(x2cs)
        X3cs.append(x3cs)
        X4cs.append(x4cs)
        ycs=sss.share(F,ycs,t,server,scale)

        mul11=shamirproc.mult2(F,Ac11,z1cs,t,server,scale)
        mul12=shamirproc.mult2(F,Ac12,z2cs,t,server,scale)
        mul13=shamirproc.mult2(F,Ac13,z3cs,t,server,scale)
        mul14=shamirproc.mult2(F,Ac14,z4cs,t,server,scale)
        mul15=shamirproc.mult2(F,Bc11,us,t,server,scale)
        mul16=shamirproc.mult2(F,Bc12,ycs,t,server,scale)

        z1cs = shamirproc.add(F,
             shamirproc.add(F,
             shamirproc.add(F,
            shamirproc.add(F,
            shamirproc.add(F,mul11,mul12,t,server,scale),mul13,t,server,scale),mul14,t,server,scale),mul15,t,server,scale),mul16,t,server,scale)

        mul21=shamirproc.mult2(F,Ac21,z1cs,t,server,scale)
        mul22=shamirproc.mult2(F,Ac22,z2cs,t,server,scale)
        mul23=shamirproc.mult2(F,Ac23,z3cs,t,server,scale)
        mul24=shamirproc.mult2(F,Ac24,z4cs,t,server,scale)
        mul25=shamirproc.mult2(F,Bc21,us,t,server,scale)
        mul26=shamirproc.mult2(F,Bc22,ycs,t,server,scale)

        z2cs = shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,mul21,mul22,t,server,scale),mul23,t,server,scale),mul24,t,server,scale),mul25,t,server,scale),mul26,t,server,scale)

        mul31=shamirproc.mult2(F,Ac31,z1cs,t,server,scale)
        mul32=shamirproc.mult2(F,Ac32,z2cs,t,server,scale)
        mul33=shamirproc.mult2(F,Ac33,z3cs,t,server,scale)
        mul34=shamirproc.mult2(F,Ac34,z4cs,t,server,scale)
        mul35=shamirproc.mult2(F,Bc31,us,t,server,scale)
        mul36=shamirproc.mult2(F,Bc32,ycs,t,server,scale)

        z3cs = shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,mul31,mul32,t,server,scale),mul33,t,server,scale),mul34,t,server,scale),mul35,t,server,scale),mul36,t,server,scale)    

        mul41=shamirproc.mult2(F,Ac41,z1cs,t,server,scale)
        mul42=shamirproc.mult2(F,Ac42,z2cs,t,server,scale)
        mul43=shamirproc.mult2(F,Ac43,z3cs,t,server,scale)
        mul44=shamirproc.mult2(F,Ac44,z4cs,t,server,scale)
        mul45=shamirproc.mult2(F,Bc41,us,t,server,scale)
        mul46=shamirproc.mult2(F,Bc42,ycs,t,server,scale)

        z4cs = shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,mul41,mul42,t,server,scale),mul43,t,server,scale),mul44,t,server,scale),mul45,t,server,scale),mul46,t,server,scale)    
        ucs=shamirproc.mult(F,Cc14,z4cs,t,server,scale)

        Z1cs.append(z1cs)
        Z2cs.append(z2cs)
        Z3cs.append(z3cs)
        Z4cs.append(z4cs)
        Ucs.append(ucs)
        Z1cs_rec.append(sss.rec(F,z1cs,rec,scale))
        Z2cs_rec.append(sss.rec(F,z2cs,rec,scale))
        Z3cs_rec.append(sss.rec(F,z3cs,rec,scale))
        Z4cs_rec.append(sss.rec(F,z4cs,rec,scale))
        Ucs_rec.append(sss.rec(F,ucs,rec,scale))
    tocs=time.time()
    print('running time of insecure controller: {0}'.format(toc-tic))
    print('running time of secure controller: {0}'.format(tocs-tics))
    print('running time of secure /running time of insecure : {0}'.format((tocs-tics)/(toc-tic)))

       ##############################################Noisy Secure controlled system######################################


    rs=sss.share(F,r,t,server,scale)
    x1csn,x2csn,x3csn,x4csn,ycsn,usn,usn_rec=0,0,0,0,0,0,0
    X1csn=[]
    X2csn=[]
    X3csn=[]
    X4csn=[]
    Ycsn=[]
    Usn_rec=[]
    Usn=[]
    ucsn=sss.share(F,0,t,server,scale)
    z1csn=sss.share(F,0,t,server,scale)
    z2csn=sss.share(F,0,t,server,scale)
    z3csn=sss.share(F,0,t,server,scale)
    z4csn=sss.share(F,0,t,server,scale)

    Z1csn=[]
    Z2csn=[]
    Z3csn=[]
    Z4csn=[]
    Ucsn=[]
    Z1csn_rec=[]
    Z2csn_rec=[]
    Z3csn_rec=[]
    Z4csn_rec=[]
    Ucsn_rec=[]

    #sss.R2F(F,5.2323,scale)

    Ac11n=sss.R2F(F,1,scale)
    Ac12n=sss.R2F(F,-1.27e-05,scale)
    Ac13n=sss.R2F(F,-0.0009439,scale)
    Ac14n=sss.R2F(F,-0.02316,scale)

    Ac21n=sss.R2F(F,0.04,scale)
    Ac22n=sss.R2F(F,1,scale)
    Ac23n=sss.R2F(F,-0.002388,scale)
    Ac24n=sss.R2F(F,-0.05874,scale)

    Ac31n=sss.R2F(F,0.0007996,scale)
    Ac32n=sss.R2F(F,+0.03996,scale)
    Ac33n=sss.R2F(F,0.997,scale)
    Ac34n=sss.R2F(F,-0.0733,scale)

    Ac41n=sss.R2F(F,2.075e-05,scale)
    Ac42n=sss.R2F(F,0.001542,scale)
    Ac43n=sss.R2F(F,0.07567,scale)
    Ac44n=sss.R2F(F,0.8928,scale)

    Bc11n=sss.R2F(F,+0.0001303,scale)    
    Bc12n=sss.R2F(F,-0.006641,scale)

    Bc21n=sss.R2F(F,+0.0001813,scale)
    Bc22n=sss.R2F(F,-0.01312,scale)

    Bc31n=sss.R2F(F,+0.0001282,scale)
    Bc32n=sss.R2F(F,-0.003186,scale)

    Bc41n=sss.R2F(F,+9.078e-05,scale)
    Bc42n=sss.R2F(F,+0.5447,scale)

    Cc14n=sss.R2F(F,64,scale)

    #Ac11n=sss.share(F,1,t,server,scale)
    #Ac12n=sss.share(F,-1.27e-05,t,server,scale)
    #Ac13n=sss.share(F,-0.0009439,t,server,scale)
    #Ac14n=sss.share(F,-0.02316,t,server,scale)

    #Ac21n=sss.share(F,0.04,t,server,scale)
    #Ac22n=sss.share(F,1,t,server,scale)
    #Ac23n=sss.share(F,-0.002388,t,server,scale)
    #Ac24n=sss.share(F,-0.05874,t,server,scale)

    #Ac31n=sss.share(F,0.0007996,t,server,scale)
    #Ac32n=sss.share(F,+0.03996,t,server,scale)
    #Ac33n=sss.share(F,0.997,t,server,scale)
    #Ac34n=sss.share(F,-0.0733,t,server,scale)

    #Ac41n=sss.share(F,2.075e-05,t,server,scale)
    #Ac42n=sss.share(F,0.001542,t,server,scale)
    #Ac43n=sss.share(F,0.07567,t,server,scale)
    #Ac44n=sss.share(F,0.8928,t,server,scale)

    #Bc11n=sss.share(F,-0.0001303,t,server,scale)   
    #Bc12n=sss.share(F,+0.006641,t,server,scale)

    #Bc21n=sss.share(F,-0.0001813,t,server,scale)
    #Bc22n=sss.share(F,+0.01312,t,server,scale)

    #Bc31n=sss.share(F,-0.0001282,t,server,scale)
    #Bc32n=sss.share(F,+0.003186,t,server,scale)

    #Bc41n=sss.share(F,-9.078e-05,t,server,scale)
    #Bc42n=sss.share(F,-0.5447,t,server,scale)

    #Cc14n=sss.share(F,64,t,server,scale)



    krs=sss.share(F,0.3933,t,server,scale)
    rec=sss.basispoly(F,server)

    for i in range(samp):
        usn=shamirproc.add(F,shamirproc.mult(F,rs,krs,t,server,scale),-ucsn,t,server,scale)
        Usn.append(usn)
        usn_rec=sss.rec(F,usn,rec,scale)
        Usn_rec.append(usn_rec)
        x1csn=0.9998*x1csn-0.004999*x2csn+0.004999*usn_rec
        x2csn=0.004999*x1csn+1*x2csn+1.25e-05*usn_rec
        x3csn=1.25e-05*x1csn+0.005*x2csn+1*x3csn+2.083e-08*usn_rec
        x4csn=2.083e-08*x1csn+1.25e-05*x2csn+0.005*x3csn+1*x4csn+2.604e-11*usn_rec
        ycsn=0.036*x3csn+0.9*x4csn+noise[i]
        Ycsn.append(ycsn)
        X1csn.append(x1csn)
        X2csn.append(x2csn)
        X3csn.append(x3csn)
        X4csn.append(x4csn)
        ycsn=sss.share(F,ycsn,t,server,scale)

        mul11n=shamirproc.mult2(F,Ac11n,z1csn,t,server,scale)
        mul12n=shamirproc.mult2(F,Ac12n,z2csn,t,server,scale)
        mul13n=shamirproc.mult2(F,Ac13n,z3csn,t,server,scale)
        mul14n=shamirproc.mult2(F,Ac14n,z4csn,t,server,scale)
        mul15n=shamirproc.mult2(F,Bc11n,usn,t,server,scale)
        mul16n=shamirproc.mult2(F,Bc12n,ycsn,t,server,scale)

        z1csn = shamirproc.add(F,
             shamirproc.add(F,
             shamirproc.add(F,
            shamirproc.add(F,
            shamirproc.add(F,mul11n,mul12n,t,server,scale),mul13n,t,server,scale),mul14n,t,server,scale),mul15n,t,server,scale),mul16n,t,server,scale)

        mul21n=shamirproc.mult2(F,Ac21n,z1csn,t,server,scale)
        mul22n=shamirproc.mult2(F,Ac22n,z2csn,t,server,scale)
        mul23n=shamirproc.mult2(F,Ac23n,z3csn,t,server,scale)
        mul24n=shamirproc.mult2(F,Ac24n,z4csn,t,server,scale)
        mul25n=shamirproc.mult2(F,Bc21n,usn,t,server,scale)
        mul26n=shamirproc.mult2(F,Bc22n,ycsn,t,server,scale)

        z2csn = shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,mul21n,mul22n,t,server,scale),mul23n,t,server,scale),mul24n,t,server,scale),mul25n,t,server,scale),mul26n,t,server,scale)

        mul31n=shamirproc.mult2(F,Ac31n,z1csn,t,server,scale)
        mul32n=shamirproc.mult2(F,Ac32n,z2csn,t,server,scale)
        mul33n=shamirproc.mult2(F,Ac33n,z3csn,t,server,scale)
        mul34n=shamirproc.mult2(F,Ac34n,z4csn,t,server,scale)
        mul35n=shamirproc.mult2(F,Bc31n,usn,t,server,scale)
        mul36n=shamirproc.mult2(F,Bc32n,ycsn,t,server,scale)

        z3csn = shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,mul31n,mul32n,t,server,scale),mul33n,t,server,scale),mul34n,t,server,scale),mul35n,t,server,scale),mul36n,t,server,scale)    

        mul41n=shamirproc.mult2(F,Ac41n,z1csn,t,server,scale)
        mul42n=shamirproc.mult2(F,Ac42n,z2csn,t,server,scale)
        mul43n=shamirproc.mult2(F,Ac43n,z3csn,t,server,scale)
        mul44n=shamirproc.mult2(F,Ac44n,z4csn,t,server,scale)
        mul45n=shamirproc.mult2(F,Bc41n,usn,t,server,scale)
        mul46n=shamirproc.mult2(F,Bc42n,ycsn,t,server,scale)

        z4csn = shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,
               shamirproc.add(F,mul41n,mul42n,t,server,scale),mul43n,t,server,scale),mul44n,t,server,scale),mul45n,t,server,scale),mul46n,t,server,scale)    
        ucsn=shamirproc.mult(F,Cc14n,z4csn,t,server,scale)

        Z1csn.append(z1csn)
        Z2csn.append(z2csn)
        Z3csn.append(z3csn)
        Z4csn.append(z4csn)
        Ucsn.append(ucsn)
        Z1csn_rec.append(sss.rec(F,z1csn,rec,scale))
        Z2csn_rec.append(sss.rec(F,z2csn,rec,scale))
        Z3csn_rec.append(sss.rec(F,z3csn,rec,scale))
        Z4csn_rec.append(sss.rec(F,z4csn,rec,scale))
        Ucsn_rec.append(sss.rec(F,ucsn,rec,scale))

        
        
        
        
    Us=np.array(Us,dtype='float128')
    Ucs=np.array(Ucs,dtype='float128')
    Z1cs=np.array(Z1cs,dtype='float128')
    Z2cs=np.array(Z2cs,dtype='float128')
    Z3cs=np.array(Z3cs,dtype='float128')
    Z4cs=np.array(Z4cs,dtype='float128')

    Usn=np.array(Usn,dtype='float128')
    Ucsn=np.array(Ucsn,dtype='float128')
    Z1csn=np.array(Z1csn,dtype='float128')
    Z2csn=np.array(Z2csn,dtype='float128')
    Z3csn=np.array(Z3csn,dtype='float128')
    Z4csn=np.array(Z4csn,dtype='float128')  
    

    
    control.X1c=X1c
    control.X2c=X2c
    control.X3c=X3c
    control.X4c=X4c
    
    control.X1cs=X1cs
    control.X2cs=X2cs
    control.X3cs=X3cs
    control.X4cs=X4cs
    
    control.X1csn=X1csn
    control.X2csn=X2csn
    control.X3csn=X3csn
    control.X4csn=X4csn  

    control.Z1c=Z1c
    control.Z2c=Z2c
    control.Z3c=Z3c
    control.Z4c=Z4c
    
    control.Z1cs=Z1cs
    control.Z2cs=Z2cs
    control.Z3cs=Z3cs
    control.Z4cs=Z4cs
    
    control.Z1csn=Z1csn
    control.Z2csn=Z2csn
    control.Z3csn=Z3csn
    control.Z4csn=Z4csn
    
    control.Z1cs_rec=Z1cs_rec
    control.Z2cs_rec=Z2cs_rec
    control.Z3cs_rec=Z3cs_rec
    control.Z4cs_rec=Z4cs_rec
    
    control.Z1csn_rec=Z1csn_rec
    control.Z2csn_rec=Z2csn_rec
    control.Z3csn_rec=Z3csn_rec
    control.Z4csn_rec=Z4csn_rec
    
    control.U=U
    control.Us=Us
    control.Usn=Usn
    
    control.Uc=Uc
    control.Ucs=Ucs
    control.Ucsn=Ucsn

    control.Us_rec=Us_rec
    control.Usn_rec=Usn_rec
    control.Ucs_rec=Ucs_rec
    control.Ucsn_rec=Ucsn_rec
    
    control.noise=noise
    control.samples=samples
    #signal energy and power controlled-inusEcure

    x1cE=0
    for i in abs(np.array(X1c)):
        x1cE+=i**2
    X1cP='%.3E' % Decimal(str(x1cE/samp))

    x2cE=0
    for i in abs(np.array(X2c)):
        x2cE+=i**2
    X2cP='%.3E' % Decimal(str(x2cE/samp))

    x3cE=0
    for i in abs(np.array(X3c)):
        x3cE+=i**2
    X3cP='%.3E' % Decimal(str(x3cE/samp))

    x4cE=0
    for i in abs(np.array(X4c)):
        x4cE+=i**2
    X4cP='%.3E' % Decimal(str(x4cE/samp))

    z1cE=0
    for i in abs(np.array(Z1c)):
        z1cE+=i**2
    Z1cP='%.3E' % Decimal(str(z1cE/samp))

    z2cE=0
    for i in abs(np.array(Z2c)):
        z2cE+=i**2
    Z2cP='%.3E' % Decimal(str(z2cE/samp))

    z3cE=0
    for i in abs(np.array(Z3c)):
        z3cE+=i**2
    Z3cP='%.3E' % Decimal(str(z3cE/samp))

    z4cE=0
    for i in abs(np.array(Z4c)):
        z4cE+=i**2
    Z4cP='%.3E' % Decimal(str(z4cE/samp))    

    uE=0
    for i in abs(np.array(U)):
        uE+=i**2
    uP='%.3E' % Decimal(str(uE/samp))

    ycE=0
    for i in abs(np.array(Yc)):
        ycE+=i**2
    YcP='%.3E' % Decimal(str(ycE/samp))


    #signal energy and power controlled-inusEcure noisy

    x1cnE=0
    for i in abs(np.array(X1cn)):
        x1cnE+=i**2
    X1cnP='%.3E' % Decimal(str(x1cnE/samp))

    x2cnE=0
    for i in abs(np.array(X2cn)):
        x2cnE+=i**2
    X2cnP='%.3E' % Decimal(str(x2cnE/samp))

    x3cnE=0
    for i in abs(np.array(X3cn)):
        x3cnE+=i**2
    X3cnP='%.3E' % Decimal(str(x3cnE/samp))

    x4cnE=0
    for i in abs(np.array(X4cn)):
        x4cnE+=i**2
    X4cnP='%.3E' % Decimal(str(x4cnE/samp))

    z1cnE=0
    for i in abs(np.array(Z1cn)):
        z1cnE+=i**2
    Z1cnP='%.3E' % Decimal(str(z1cnE/samp))

    z2cnE=0
    for i in abs(np.array(Z2cn)):
        z2cnE+=i**2
    Z2cnP='%.3E' % Decimal(str(z2cnE/samp))

    z3cnE=0
    for i in abs(np.array(Z3cn)):
        z3cnE+=i**2
    Z3cnP='%.3E' % Decimal(str(z3cnE/samp))

    z4cnE=0
    for i in abs(np.array(Z4cn)):
        z4cnE+=i**2
    Z4cnP='%.3E' % Decimal(str(z4cnE/samp))    

    unE=0
    for i in abs(np.array(Un)):
        unE+=i**2
    unP='%.3E' % Decimal(str(unE/samp))

    ycnE=0
    for i in abs(np.array(Ycn)):
        ycnE+=i**2
    YcnP='%.3E' % Decimal(str(ycnE/samp))

    #signal energy and power controlled-sEcure

    X1cusE=0
    for i in abs(np.array(X1cs)):
        X1cusE+=i**2
    X1csP='%.3E' % Decimal(str(X1cusE/samp))

    X2cusE=0
    for i in abs(np.array(X2cs)):
        X2cusE+=i**2
    X2csP='%.3E' % Decimal(str(X2cusE/samp))

    X3cusE=0
    for i in abs(np.array(X3cs)):
        X3cusE+=i**2
    X3csP='%.3E' % Decimal(str(X3cusE/samp))

    X4cusE=0
    for i in abs(np.array(X4cs)):
        X4cusE+=i**2
    X4csP='%.3E' % Decimal(str(X4cusE/samp))    

    z1csE=0
    for i in abs(np.array(Z1cs_rec)):
        z1csE+=i**2
    Z1csP='%.3E' % Decimal(str(z1csE/samp))

    z2csE=0
    for i in abs(np.array(Z2cs_rec)):
        z2csE+=i**2
    Z2csP='%.3E' % Decimal(str(z2csE/samp))

    z3csE=0
    for i in abs(np.array(Z3cs_rec)):
        z3csE+=i**2
    Z3csP='%.3E' % Decimal(str(z3csE/samp))

    z4csE=0
    for i in abs(np.array(Z4cs_rec)):
        z4csE+=i**2
    Z4csP='%.3E' % Decimal(str(z4csE/samp))    

    usE=0
    for i in abs(np.array(Us_rec)):
        usE+=i**2
    usP='%.3E' % Decimal(str(usE/samp))

    ycsE=0
    for i in abs(np.array(Ycs)):
        ycsE+=i**2
    YcsP='%.3E' % Decimal(str(ycsE/samp))


    #signal energy and power noisy controlled-sEcure

    X1cusnE=0
    for i in abs(np.array(X1csn)):
        X1cusnE+=i**2
    X1csnP='%.3E' % Decimal(str(X1cusnE/samp))

    X2cusnE=0
    for i in abs(np.array(X2csn)):
        X2cusnE+=i**2
    X2csnP='%.3E' % Decimal(str(X2cusnE/samp))

    X3cusnE=0
    for i in abs(np.array(X3csn)):
        X3cusnE+=i**2
    X3csnP='%.3E' % Decimal(str(X3cusnE/samp))

    X4cusnE=0
    for i in abs(np.array(X4csn)):
        X4cusnE+=i**2
    X4csnP='%.3E' % Decimal(str(X4cusnE/samp))    

    z1csnE=0
    for i in abs(np.array(Z1csn_rec)):
        z1csnE+=i**2
    Z1csnP='%.3E' % Decimal(str(z1csnE/samp))

    z2csnE=0
    for i in abs(np.array(Z2csn_rec)):
        z2csnE+=i**2
    Z2csnP='%.3E' % Decimal(str(z2csnE/samp))

    z3csnE=0
    for i in abs(np.array(Z3csn_rec)):
        z3csnE+=i**2
    Z3csnP='%.3E' % Decimal(str(z3csnE/samp))

    z4csnE=0
    for i in abs(np.array(Z4csn_rec)):
        z4csnE+=i**2
    Z4csnP='%.3E' % Decimal(str(z4csnE/samp))    

    usnE=0
    for i in abs(np.array(Usn_rec)):
        usnE+=i**2
    usnP='%.3E' % Decimal(str(usnE/samp))

    ycsnE=0
    for i in abs(np.array(Ycsn)):
        ycsnE+=i**2
    YcsnP='%.3E' % Decimal(str(ycsnE/samp))    


    ##differences
    x1diffrence = abs(np.array(X1c)-np.array(X1cs))
    x2diffrence = abs(np.array(X2c)-np.array(X2cs))
    x3diffrence = abs(np.array(X3c)-np.array(X3cs))
    x4diffrence = abs(np.array(X4c)-np.array(X4cs))

    z1diffrence = abs(np.array(Z1c)-np.array(Z1cs_rec))
    z2diffrence = abs(np.array(Z2c)-np.array(Z2cs_rec))
    z3diffrence = abs(np.array(Z3c)-np.array(Z3cs_rec))
    z4diffrence = abs(np.array(Z4c)-np.array(Z4cs_rec))
    ydiffrence = abs(np.array(Yc)-np.array(Ycs))
    udiffrence = abs(np.array(U)-np.array(Us_rec))

    ##differences noisy
    x1ndiffrence = abs(np.array(X1cn)-np.array(X1csn))
    x2ndiffrence = abs(np.array(X2cn)-np.array(X2csn))
    x3ndiffrence = abs(np.array(X3cn)-np.array(X3csn))
    x4ndiffrence = abs(np.array(X4cn)-np.array(X4csn))

    z1ndiffrence = abs(np.array(Z1cn)-np.array(Z1csn_rec))
    z2ndiffrence = abs(np.array(Z2cn)-np.array(Z2csn_rec))
    z3ndiffrence = abs(np.array(Z3cn)-np.array(Z3csn_rec))
    z4ndiffrence = abs(np.array(Z4cn)-np.array(Z4csn_rec))
    yndiffrence = abs(np.array(Ycn)-np.array(Ycsn))
    undiffrence = abs(np.array(Un)-np.array(Usn_rec))

    #table of differences (without noise)
    sx1=0
    for i in x1diffrence:
        sx1+=i**2
    IAEx1= '%.3E' % Decimal(str(sum(x1diffrence)/samp))
    ISEx1= '%.3E' % Decimal(str(sx1/samp))

    sx2=0
    for i in x2diffrence:
        sx2+=i**2
    IAEx2= '%.3E' % Decimal(str(sum(x2diffrence)/samp))
    ISEx2= '%.3E' % Decimal(str(sx2/samp))

    sx3=0
    for i in x3diffrence:
        sx3+=i**2
    IAEx3= '%.3E' % Decimal(str(sum(x3diffrence)/samp))
    ISEx3= '%.3E' % Decimal(str(sx3/samp))

    sx4=0
    for i in x4diffrence:
        sx4+=i**2
    IAEx4= '%.3E' % Decimal(str(sum(x4diffrence)/samp))
    ISEx4= '%.3E' % Decimal(str(sx4/samp))

    sz1=0
    for i in z1diffrence:
        sz1+=i**2
    IAEz1= '%.3E' % Decimal(str(sum(z1diffrence)/samp))
    ISEz1= '%.3E' % Decimal(str(sz1/samp))

    sz2=0
    for i in z2diffrence:
        sz2+=i**2
    IAEz2= '%.3E' % Decimal(str(sum(z2diffrence)/samp))
    ISEz2= '%.3E' % Decimal(str(sz2/samp))   

    sz3=0
    for i in z3diffrence:
        sz3+=i**2
    IAEz3= '%.3E' % Decimal(str(sum(z3diffrence)/samp))
    ISEz3= '%.3E' % Decimal(str(sz3/samp))    

    sz4=0
    for i in z1diffrence:
        sz4+=i**2
    IAEz4= '%.3E' % Decimal(str(sum(z4diffrence)/samp))
    ISEz4= '%.3E' % Decimal(str(sz4/samp))

    su=0
    for i in udiffrence:
        su+=i**2
    IAEsu= '%.3E' % Decimal(str(sum(udiffrence)/samp))
    ISEsu= '%.3E' % Decimal(str(su/samp))

    sy=0
    for i in ydiffrence:
        sy+=i**2
    IAEsy= '%.3E' % Decimal(str(sum(ydiffrence)/samp))
    ISEsy= '%.3E' % Decimal(str(sy/samp))



    table = [['Signal (without noise)','IAE','ISE','CL/Inusecure','CL/Secure']
             ,["x1",IAEx1,ISEx1,X1cP,X1csP]
            ,["x2",IAEx2,ISEx2,X2cP,X2csP]
             ,["x3",IAEx3,ISEx3,X3cP,X3csP]
          ,["x4",IAEx4,ISEx4,X4cP,X4csP]
             ,["z1",IAEz1,ISEz1,Z1cP,Z1csP]
             ,["z2",IAEz2,ISEz2,Z2cP,Z2csP]
            ,["z3",IAEz3,ISEz3,Z3cP,Z3csP]
            ,["z4",IAEz4,ISEz4,Z4cP,Z4csP]
            ,["u",IAEsu,ISEsu,uP,usP]
            ,["y (theta2)",IAEsy,ISEsy,YcP,YcsP]]

    headers = ['Signal (without noise)',"Integral Absolute Error(IAE)\nsum[abs(inusEcuretheta-usEcuretheta)]/samples", "Integral Squared Error(IusE)\nsum[abs(sqrt(inusEcuretheta-usEcuretheta))]/samples",'ClousEd Loop/inusEcure','ClousEd Loop/usEcure']

    print(tabulate(table, headers, tablefmt="pretty", disable_numparse=True))
    fold=open('Results(without noise).txt','w')
    fold.write(tabulate(table))
    fold.close()


    #table of differences (with noise)
    sx1n=0
    for i in x1ndiffrence:
        sx1n+=i**2
    IAEx1n= '%.3E' % Decimal(str(sum(x1ndiffrence)/samp))
    ISEx1n= '%.3E' % Decimal(str(sx1n/samp))

    sx2n=0
    for i in x2ndiffrence:
        sx2n+=i**2
    IAEx2n= '%.3E' % Decimal(str(sum(x2ndiffrence)/samp))
    ISEx2n= '%.3E' % Decimal(str(sx2n/samp))

    sx3n=0
    for i in x3ndiffrence:
        sx3n+=i**2
    IAEx3n= '%.3E' % Decimal(str(sum(x3ndiffrence)/samp))
    ISEx3n= '%.3E' % Decimal(str(sx3n/samp))

    sx4n=0
    for i in x4ndiffrence:
        sx4n+=i**2
    IAEx4n= '%.3E' % Decimal(str(sum(x4ndiffrence)/samp))
    ISEx4n= '%.3E' % Decimal(str(sx4n/samp))

    sz1n=0
    for i in z1ndiffrence:
        sz1n+=i**2
    IAEz1n= '%.3E' % Decimal(str(sum(z1ndiffrence)/samp))
    ISEz1n= '%.3E' % Decimal(str(sz1n/samp))

    sz2n=0
    for i in z2ndiffrence:
        sz2n+=i**2
    IAEz2n= '%.3E' % Decimal(str(sum(z2ndiffrence)/samp))
    ISEz2n= '%.3E' % Decimal(str(sz2n/samp))   

    sz3n=0
    for i in z3ndiffrence:
        sz3n+=i**2
    IAEz3n= '%.3E' % Decimal(str(sum(z3ndiffrence)/samp))
    ISEz3n= '%.3E' % Decimal(str(sz3n/samp))    

    sz4n=0
    for i in z1ndiffrence:
        sz4n+=i**2
    IAEz4n= '%.3E' % Decimal(str(sum(z4ndiffrence)/samp))
    ISEz4n= '%.3E' % Decimal(str(sz4n/samp))

    sun=0
    for i in undiffrence:
        sun+=i**2
    IAEsun= '%.3E' % Decimal(str(sum(undiffrence)/samp))
    ISEsun= '%.3E' % Decimal(str(sun/samp))

    syn=0
    for i in yndiffrence:
        syn+=i**2
    IAEsyn= '%.3E' % Decimal(str(sum(yndiffrence)/samp))
    ISEsyn= '%.3E' % Decimal(str(syn/samp))



    table2 = [['Signal (with noise)','IAE','ISE','CL/Inusecure','CL/Secure']
             ,["x1",IAEx1n,ISEx1n,X1cnP,X1csnP]
            ,["x2",IAEx2n,ISEx2n,X2cnP,X2csnP]
             ,["x3",IAEx3n,ISEx3n,X3cnP,X3csnP]
          ,["x4",IAEx4n,ISEx4n,X4cnP,X4csnP]
             ,["z1",IAEz1n,ISEz1n,Z1cnP,Z1csnP]
             ,["z2",IAEz2n,ISEz2n,Z2cnP,Z2csnP]
            ,["z3",IAEz3n,ISEz3n,Z3cnP,Z3csnP]
            ,["z4",IAEz4n,ISEz4n,Z4cnP,Z4csnP]
            ,["u",IAEsun,ISEsun,unP,usnP]
            ,["y (theta2)",IAEsyn,ISEsyn,YcnP,YcsnP]]

    headers2 = ['Signal (with noise)',"Integral Absolute Error(IAE)\nsum[abs(inusEcuretheta-usEcuretheta)]/samples", "Integral Squared Error(IusE)\nsum[abs(sqrt(inusEcuretheta-usEcuretheta))]/samples",'ClousEd Loop/inusEcure','ClousEd Loop/usEcure']

    print(tabulate(table2, headers2, tablefmt="pretty", disable_numparse=True))
    fold2=open('Results(with noise).txt','w')
    fold2.write(tabulate(table2))
    fold2.close()    




    ##INSECURE CONTROLLED SYSTEM STEP RESPONSE CHARACTERISTIC##
    print('\n\033[1m'+'INSECURE CONTROLLED SYSTEM STEP RESPONSE CHARACTERISTICS\n'+'\033[0m')
    fcont = interp1d(Yc,samples) 
    f2cont=interp1d(samples,Yc) 
    #stepresponse of sys
    contau = fcont(0.632*r)
    contrising = fcont(0.85*r)-fcont(0.05*r)
    contsettling = fcont(0.98*r)
    tau=print('Insecure controlled system constant time :{0}'.format(contau))
    risingtime = print('Insecure controlled system rise_time :{0}'.format(contrising)) 
    settlingtime = print('Insecure controlled system settlingtime :{0}'.format(contsettling))
    if max(Yc)>r:
        peakamp=round(max(Yc),2)
        overshoot=round(((max(Yc)-r)/r)*100,1)
        peaktime=round(float(fcont(max(Yc))),2)
        print('Insecure controlled system "peak amplitude :{0}" "overshoot :{1}" "At time :{2}"'
              .format(peakamp,overshoot,peaktime))

    ##INSECURE NOISY CONTROLLED SYSTEM STEP RESPONSE CHARACTERISTIC##
    print('\n\033[1m'+'INSECURE NOISY CONTROLLED SYSTEM STEP RESPONSE CHARACTERISTICS\n'+'\033[0m')
    fncont = interp1d(Ycn,samples) 
    f2ncont=interp1d(samples,Ycn) 
    #stepresponse of sys
    contaun = fncont(0.632*r)
    contrisingn = fncont(0.85*r)-fncont(0.05*r)
    contsettlingn = fncont(0.98*r)
    tau=print('Insecure noisy controlled system constant time :{0}'.format(contaun))
    risingtime = print('Insecure noisy controlled system rise_time :{0}'.format(contrisingn)) 
    settlingtime = print('Insecure noisy controlled system settlingtime :{0}'.format(contsettlingn))
    if max(Ycn)>r:
        peakampn=round(max(Ycn),2)
        overshootn=round(((max(Ycn)-r)/r)*100,1)
        peaktimen=round(float(fncont(max(Ycn))),2)
        print('Insecure noisy controlled system "peak amplitude :{0}" "overshoot :{1}" "At time :{2}"'
              .format(peakampn,overshootn,peaktimen))    

    ##SECURE STEP CONTROLLED SYSTEM RESPONSE CHARACTERISTIC##
    print('\n\033[1m'+'SECURE CONTROLLED SYSTEM STEP RESPONSE CHARACTERISTICS\n'+'\033[0m')
    fscont = interp1d(Ycs,samples)
    f2scont = interp1d(samples,Ycs)
    #stepresponse of sys
    scontau = fscont(0.632*r)
    scontrising = fscont(0.85*r)-fscont(0.05*r)
    scontsettling = fscont(0.98*r)
    tau=print('Secure controlled system constant time :{0}'.format(scontau))
    risingtime = print('Secure controlled system rise_time :{0}'.format(scontrising)) 
    settlingtime = print('Secure controlled system settlingtime :{0}'.format(scontsettling))
    if max(Ycs)>r:
        peakamps=round(max(Ycs)/r,2)
        overshoots=round(((max(Ycs)-r)/r)*100,1)
        peaktimes=round(float(fscont(max(Ycs))),2)
        print('Secure controlled system "peak amplitude :{0}" "overshoot :{1}" "At time :{2}"'
              .format(peakamps,overshoots,peaktimes))    

    ##SECURE NOISY STEP CONTROLLED SYSTEM RESPONSE CHARACTERISTIC##
    print('\n\033[1m'+'SECURE NOISY CONTROLLED SYSTEM STEP RESPONSE CHARACTERISTICS\n'+'\033[0m')
    fsncont = interp1d(Ycsn,samples)
    f2sncont = interp1d(samples,Ycsn)
    #stepresponse of sys
    scontaun = fsncont(0.632*r)
    scontrisingn = fsncont(0.85*r)-fsncont(0.05*r)
    scontsettlingn = fsncont(0.98*r)
    tau=print('Secure noisy controlled system constant time :{0}'.format(scontaun))
    risingtime = print('Secure noisy controlled system rise_time :{0}'.format(scontrisingn)) 
    settlingtime = print('Secure noisy controlled system settlingtime :{0}'.format(scontsettlingn))
    if max(Ycsn)>r:
        peakampsn=round(max(Ycsn)/r,2)
        overshootsn=round(((max(Ycsn)-r)/r)*100,1)
        peaktimesn=round(float(fsncont(max(Ycsn))),2)
        print('Secure noisy controlled system "peak amplitude :{0}" "overshoot :{1}" "At time :{2}"'
              .format(peakampsn,overshootsn,peaktimesn)) 

    ####################savemat#######################    
    Us=np.array(Us,dtype='float128')
    Ucs=np.array(Ucs,dtype='float128')
    Z1cs=np.array(Z1cs,dtype='float128')
    Z2cs=np.array(Z2cs,dtype='float128')
    Z3cs=np.array(Z3cs,dtype='float128')
    Z4cs=np.array(Z4cs,dtype='float128')

    Usn=np.array(Usn,dtype='float128')
    Ucsn=np.array(Ucsn,dtype='float128')
    Z1csn=np.array(Z1csn,dtype='float128')
    Z2csn=np.array(Z2csn,dtype='float128')
    Z3csn=np.array(Z3csn,dtype='float128')
    Z4csn=np.array(Z4csn,dtype='float128')

    ###save matlab
    results={
        'X1':X1,
        'X2':X2,
        'X3':X3,
        'X4':X4,
        'Y':Y,

        'X1c':X1c,
        'X2c':X2c,
        'X3c':X3c,
        'X4c':X4c,
        'Yc':Yc,
        'Z1c':Z1c,
        'Z2c':Z2c,
        'Z3c':Z3c,
        'Z4c':Z4c,
        'U':U,
        'Uc':Uc,
        
        'X1cn':X1cn,
        'X2cn':X2cn,
        'X3cn':X3cn,
        'X4cn':X4cn,
        'Ycn':Ycn,
        'Z1cn':Z1cn,
        'Z2cn':Z2cn,
        'Z3cn':Z3cn,
        'Z4cn':Z4cn,
        'Un':Un,
        'Ucn':Ucn,
        
        'X1cs':X1cs,
        'X2cs':X2cs,
        'X3cs':X3cs,
        'X4cs':X4cs,
        'Ycs':Ycs,
        'Us':Us,
        'Ucs':Ucs,
        'Us_rec':Us_rec,
        'Ucs_rec':Ucs_rec,    
        'Z1cs':Z1cs,
        'Z2cs':Z2cs,
        'Z3cs':Z3cs,
        'Z4cs':Z4cs,
        'Z1cs_rec':Z1cs_rec,
        'Z2cs_rec':Z2cs_rec,
        'Z3cs_rec':Z3cs_rec,
        'Z4cs_rec':Z4cs_rec,

        'X1csn':X1csn,
        'X2csn':X2csn,
        'X3csn':X3csn,
        'X4csn':X4csn,
        'Ycsn':Ycsn,
        'Usn':Usn,
        'Ucsn':Ucsn,
        'Usn_rec':Usn_rec,
        'Ucsn_rec':Ucsn_rec,    
        'Z1csn':Z1csn,
        'Z2csn':Z2csn,
        'Z3csn':Z3csn,
        'Z4csn':Z4csn,
        'Z1csn_rec':Z1csn_rec,
        'Z2csn_rec':Z2csn_rec,
        'Z3csn_rec':Z3csn_rec,
        'Z4csn_rec':Z4csn_rec,
        'noise':noise,
        'samples':samples
        }
    savemat("results.mat",results)    


    ###############################################plots########################################
    #insecure controlled sys
    fig,axs = plt.subplots(2, 2)
    fig.set_size_inches(20, 13, forward=True)
    #fig.suptitle('Controlled System')
    axs[0, 0].step(samples,Yc,color = 'red')
    if stepinfo==True:
        axs[0, 0].scatter(contrising,f2cont(contrising),color = 'blue')
        axs[0, 0].scatter(contsettling,f2cont(contsettling),color = 'blue')
        axs[0, 0].text(contrising,f2cont(contrising),
                 '  sys:Insecure controlled\n  RisingTime : {0}'
                       .format(round(float(contrising),2)),horizontalalignment='left'
                       ,verticalalignment='bottom',fontsize=11)
        axs[0, 0].text(contsettling,f2cont(contsettling),
                 '  sys:Insecure controlled\n  Settlingtime : {0}'
                       .format(round(float(contsettling),2)),verticalalignment='top',fontsize=11)
        if max(Yc)>r:
            axs[0, 0].scatter(peaktime,peakamp,color = 'blue')
            axs[0, 0].text(peaktime,peakamp,
                 '  sys:Insecure controlled\n  Peak amplitude : {0}\n  Overshoot : {1}\n  At time : {2}'
                     .format(peakamp,overshoot,peaktime),horizontalalignment='left'
                           ,verticalalignment='top',fontsize=11)
    axs[0, 0].set_title('Insecure Controlled System')
    axs[0, 0].axhline(y = r, color = 'black', linestyle = '--')
    axs[0, 0].label_outer()

    ##insecure noisy controlled
    axs[1, 0].step(samples,Ycn,color = 'red')
    if stepinfo==True:
        axs[1, 0].scatter(contrisingn,f2ncont(contrisingn),color = 'blue')
        axs[1, 0].scatter(contsettlingn,f2ncont(contsettlingn),color = 'blue')
        axs[1, 0].text(contrisingn,f2ncont(contrisingn),
                 '  sys:Insecure noisy controlled\n  RisingTime : {0}'
                       .format(round(float(contrisingn),2)),horizontalalignment='left'
                       ,verticalalignment='bottom',fontsize=11)
        axs[1, 0].text(contsettlingn,f2ncont(contsettlingn),
                 '  sys:Insecure noisy controlled\n  Settlingtime : {0}'
                       .format(round(float(contsettlingn),2)),verticalalignment='top',fontsize=11)
        if max(Ycn)>r:
            axs[1, 0].scatter(peaktimen,peakampn,color = 'blue')
            axs[1, 0].text(peaktimen,peakampn,
                 '  sys:Insecure noisy controlled\n  Peak amplitude : {0}\n  Overshoot : {1}\n  At time : {2}'
                     .format(peakampn,overshootn,peaktimen),horizontalalignment='left'
                           ,verticalalignment='top',fontsize=11)
    axs[1, 0].set_title('Insecure noisy Controlled System')
    axs[1, 0].axhline(y = r, color = 'black', linestyle = '--')

    ###secure controlled

    axs[0, 1].step(samples,Ycs,color = 'blue')
    if stepinfo==True:
        axs[0, 1].scatter(scontrising,f2scont(scontrising),color = 'red')
        axs[0, 1].scatter(scontsettling,f2scont(scontsettling),color = 'red')
        axs[0, 1].text(scontrising,f2scont(scontrising),
                 '  sys:Secure controlled\n  RisingTime : {0}'
                       .format(round(float(scontrising),2)),horizontalalignment='left'
                       ,verticalalignment='bottom',fontsize=11)
        axs[0, 1].text(scontsettling,f2scont(scontsettling),
                 '  sys:Secure controlled\n  Settlingtime : {0}'
                       .format(round(float(scontsettling),2)),verticalalignment='top',fontsize=11)
        if max(Ycs)>r:
            axs[0, 1].scatter(peaktimes,peakamps,color = 'red')
            axs[0, 1].text(peaktimes,peakamps,
                 '  sys:Secure controlled\n  Peak amplitude : {0}\n  Overshoot : {1}\n  At time : {2}'
                     .format(peakamps,overshoots,peaktimes),horizontalalignment='left'
                           ,verticalalignment='top',fontsize=11)
    axs[0,1].set_title('Secure Controlled System')
    axs[0, 1].axhline(y = r, color = 'black', linestyle = '--')

    ###secure noisy controlled

    axs[1, 1].step(samples,Ycsn,color = 'blue')
    if stepinfo==True:
        axs[1, 1].scatter(scontrisingn,f2sncont(scontrisingn),color = 'red')
        axs[1, 1].scatter(scontsettlingn,f2sncont(scontsettlingn),color = 'red')
        axs[1, 1].text(scontrisingn,f2sncont(scontrisingn),
                 '  sys:Secure noisy controlled\n  RisingTime : {0}'
                       .format(round(float(scontrisingn),2)),horizontalalignment='left'
                       ,verticalalignment='bottom',fontsize=11)
        axs[1, 1].text(scontsettlingn,f2sncont(scontsettlingn),
                 '  sys:Secure noisy controlled\n  Settlingtime : {0}'
                       .format(round(float(scontsettlingn),2)),verticalalignment='top',fontsize=11)
        if max(Ycsn)>r:
            axs[1, 1].scatter(peaktimesn,peakampsn,color = 'red')
            axs[1, 1].text(peaktimesn,peakampsn,
                 '  sys:Secure noisy controlled\n  Peak amplitude : {0}\n  Overshoot : {1}\n  At time : {2}'
                     .format(peakampsn,overshootsn,peaktimesn),horizontalalignment='left'
                           ,verticalalignment='top',fontsize=11)
    axs[1,1].set_title('Secure noisy Controlled System')
    axs[1, 1].axhline(y = r, color = 'black', linestyle = '--')


    ########## control signal #######
    fig1,axs1 = plt.subplots(4, 1)
    fig1.set_size_inches(20, 13, forward=True)
    axs1[0].step(samples,U,color = 'red')
    axs1[0].set_ylabel('insecure control signal(rad)',fontsize=16)

    axs1[1].step(samples,Us_rec,color = 'blue')
    axs1[1].set_ylabel('secure control signal(rad)',fontsize=16)
    axs1[1].set_xlabel('Time (sec)',fontsize=16)

    axs1[2].step(samples,Un,color = 'red')
    axs1[2].set_ylabel('insecure noisy control signal(rad)',fontsize=16)
    axs1[2].set_xlabel('Time (sec)',fontsize=16)

    axs1[3].step(samples,Usn_rec,color = 'blue')
    axs1[3].set_ylabel('secure noisy control signal(rad)',fontsize=16)
    axs1[3].set_xlabel('Time (sec)',fontsize=16)

    ########## insecure plant states #######
    fig2,axs2 = plt.subplots(2, 2)
    fig2.set_size_inches(20, 13, forward=True)
    axs2[0,0].step(samples,X1c,color = 'r')
    axs2[0,0].set_ylabel('insecure x1',fontsize=16)

    axs2[0,1].step(samples,X2c,color = 'r')
    axs2[0,1].set_ylabel('insecure x2',fontsize=16)
    axs2[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs2[1,0].step(samples,X3c,color = 'r')
    axs2[1,0].set_ylabel('insecure x3',fontsize=16)
    axs2[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs2[1,1].step(samples,X4c,color = 'r')
    axs2[1,1].set_ylabel('insecure x4',fontsize=16)
    axs2[1,1].set_xlabel('Time (sec)',fontsize=16)

    ########## secure plant states #######
    fig3,axs3 = plt.subplots(2, 2)
    fig3.set_size_inches(20, 13, forward=True)
    axs3[0,0].step(samples,X1cs,color = 'b')
    axs3[0,0].set_ylabel('secure x1',fontsize=16)

    axs3[0,1].step(samples,X2cs,color = 'b')
    axs3[0,1].set_ylabel('secure x2',fontsize=16)
    axs3[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs3[1,0].step(samples,X3cs,color = 'b')
    axs3[1,0].set_ylabel('secure x3',fontsize=16)
    axs3[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs3[1,1].step(samples,X4cs,color = 'b')
    axs3[1,1].set_ylabel('secure x4',fontsize=16)
    axs3[1,1].set_xlabel('Time (sec)',fontsize=16)

    ########## insecure noisy plant states #######
    fig4,axs4 = plt.subplots(2, 2)
    fig4.set_size_inches(20, 13, forward=True)
    axs4[0,0].step(samples,X1cn,color = 'r')
    axs4[0,0].set_ylabel('insecure noisy x1',fontsize=16)

    axs4[0,1].step(samples,X2cn,color = 'r')
    axs4[0,1].set_ylabel('insecure noisy x2',fontsize=16)
    axs4[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs4[1,0].step(samples,X3cn,color = 'r')
    axs4[1,0].set_ylabel('insecure noisy x3',fontsize=16)
    axs4[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs4[1,1].step(samples,X4cn,color = 'r')
    axs4[1,1].set_ylabel('insecure noisy x4',fontsize=16)
    axs4[1,1].set_xlabel('Time (sec)',fontsize=16)

    ########## noisy secure plant states #######
    fig5,axs5 = plt.subplots(2, 2)
    fig5.set_size_inches(20, 13, forward=True)
    axs5[0,0].step(samples,X1csn,color = 'b')
    axs5[0,0].set_ylabel('secure noisy x1',fontsize=16)

    axs5[0,1].step(samples,X2csn,color = 'b')
    axs5[0,1].set_ylabel('secure noisy x2',fontsize=16)
    axs5[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs5[1,0].step(samples,X3csn,color = 'b')
    axs5[1,0].set_ylabel('secure noisy x3',fontsize=16)
    axs5[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs5[1,1].step(samples,X4csn,color = 'b')
    axs5[1,1].set_ylabel('secure noisy x4',fontsize=16)
    axs5[1,1].set_xlabel('Time (sec)',fontsize=16)


    ########## insecure controller states #######
    fig6,axs6 = plt.subplots(2, 2)
    fig6.set_size_inches(20, 13, forward=True)
    axs6[0,0].step(samples,Z1c,color = 'r')
    axs6[0,0].set_ylabel('insecure z1',fontsize=16)

    axs6[0,1].step(samples,Z2c,color = 'r')
    axs6[0,1].set_ylabel('insecure z2',fontsize=16)
    axs6[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs6[1,0].step(samples,Z3c,color = 'r')
    axs6[1,0].set_ylabel('insecure z3',fontsize=16)
    axs6[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs6[1,1].step(samples,Z4c,color = 'r')
    axs6[1,1].set_ylabel('insecure z4',fontsize=16)
    axs6[1,1].set_xlabel('Time (sec)',fontsize=16)

    ########## secure controller states #######
    fig7,axs7 = plt.subplots(2, 2)
    fig7.set_size_inches(20, 13, forward=True)
    axs7[0,0].step(samples,Z1cs_rec,color = 'b')
    axs7[0,0].set_ylabel('secure z1',fontsize=16)

    axs7[0,1].step(samples,Z2cs_rec,color = 'b')
    axs7[0,1].set_ylabel('secure z2',fontsize=16)
    axs7[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs7[1,0].step(samples,Z3cs_rec,color = 'b')
    axs7[1,0].set_ylabel('secure z3',fontsize=16)
    axs7[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs7[1,1].step(samples,Z4cs_rec,color = 'b')
    axs7[1,1].set_ylabel('secure z4',fontsize=16)
    axs7[1,1].set_xlabel('Time (sec)',fontsize=16)

    ########## noisy insecure controller states #######
    fig8,axs8 = plt.subplots(2, 2)
    fig8.set_size_inches(20, 13, forward=True)
    axs8[0,0].step(samples,Z1cn,color = 'r')
    axs8[0,0].set_ylabel('insecure noisy z1',fontsize=16)

    axs8[0,1].step(samples,Z2cn,color = 'r')
    axs8[0,1].set_ylabel('insecure noisy z2',fontsize=16)
    axs8[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs8[1,0].step(samples,Z3cn,color = 'r')
    axs8[1,0].set_ylabel('insecure noisy z3',fontsize=16)
    axs8[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs8[1,1].step(samples,Z4cn,color = 'r')
    axs8[1,1].set_ylabel('insecure noisy z4',fontsize=16)
    axs8[1,1].set_xlabel('Time (sec)',fontsize=16)

    ########## noisy secure controller states #######
    fig9,axs9 = plt.subplots(2, 2)
    fig9.set_size_inches(20, 13, forward=True)
    axs9[0,0].step(samples,Z1csn_rec,color = 'b')
    axs9[0,0].set_ylabel('secure noisy z1',fontsize=16)

    axs9[0,1].step(samples,Z2csn_rec,color = 'b')
    axs9[0,1].set_ylabel('secure noisy z2',fontsize=16)
    axs9[0,1].set_xlabel('Time (sec)',fontsize=16)

    axs9[1,0].step(samples,Z3csn_rec,color = 'b')
    axs9[1,0].set_ylabel('secure noisy z3',fontsize=16)
    axs9[1,0].set_xlabel('Time (sec)',fontsize=16)

    axs9[1,1].step(samples,Z4csn_rec,color = 'b')
    axs9[1,1].set_ylabel('secure noisy z4',fontsize=16)
    axs9[1,1].set_xlabel('Time (sec)',fontsize=16)


    ###########server control signal computation###########
    fig10,axs10 = plt.subplots()
    fig10.set_size_inches(20, 13, forward=True)
    for i in range(server):
        axs10.scatter(samples,Us[:,i])
    axs10.set_xlim(0,Tf)
    axs10.set_ylabel('Distributed control force (N)')
    axs10.set_xlabel('Time (sec)')
    leg=[]
    for i in range(server):
        leg.append('Server {0}'.format(str(i+1)))
    axs10.legend(leg, loc='upper right')