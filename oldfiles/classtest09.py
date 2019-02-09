#same as classtest04.py, but now adds a child class with open BC

import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from matplotlib import animation

class Oscillator:

    def __init__(self,
                 numberMasses = 5,
                 relAmp = [1.0 , 0.8 , 0.6 , 0.4 , 0.4],
                 howLong = 4.0*np.pi*5.0,
                 howMany = 4000):
        self.numberMasses = numberMasses
        self.relAmp = np.array(relAmp)
        self.howLong = howLong
        self.howMany = howMany

    def setTime(self):
        time = np.linspace(0 , self.howLong , self.howMany)
        return time

    def setMatrix(self,BC):
        N = self.numberMasses
        Dntemp = []
        for ii in range(N):
            Dntemp.append(np.zeros(N))
            Dn = np.matrix(Dntemp)
        for ii in range(N):
            Dn[ii,ii] = 2.0
        for ii in range(N-1):
            Dn[ii+1,ii] = -1.0
            Dn[ii,ii+1] = -1.0
        if BC=="open":
            Dn[0,0] = 1.0
            Dn[N-1,N-1] = 1.0
        return Dn





    def diagonalizeMatrix(self,Dn,BC):
        #print(Dn)
        N = self.numberMasses  #control should go to matrixdimensions, not the number of masses, coz open BC has N-1 dimensional matrix
        eValues,eVectors=LA.eig(Dn)
        idx = np.argsort(eValues)
        eValues = eValues[idx]
        eVectors = eVectors[:,idx]
        #scale the frequencies
        if BC=="fixed":
            eValuesMin = eValues[0]
        if BC=="open":
            eValuesMin = eValues[1]
        frequencies = []
        for ii in range(N):
            frequencies.append(eValues[ii]/eValuesMin)
        # extract the amplitudes:
        amplitudes = []
        for ii in range(N):
            Atemp = []
            for jj in range(N):
                Atemp.append(self.relAmp[ii]*eVectors[jj,ii])
            amplitudes.append(Atemp)
        return frequencies , amplitudes
    
    

    def find_x_coord(self,frequencies,amplitudes,x0):
        N = self.numberMasses
        tabLength = float(N)+1.0              # length of the table
        springLength = tabLength/(float(N+1)) # unstretched length of spring
        #x0 = x_eqb
        #for ii in range(N):
        #    x0.append((-tabLength*0.5)+(ii+1)*springLength) #equilibrium positions of the masses

        #frequencies , amplitudes  = self.diagonalizeMatrix()
        #calculate the x-coordinates
        time = self.setTime()    #np.linspace(0 , self.howLong , self.howMany)
        xall = []
        for ii in range(N):
            for jj in range(N):
                xall.append(amplitudes[ii][jj]*np.cos(frequencies[ii]*time) + x0[jj])
        #time = np.linspace(0 , self.howLong , self.howMany)
        return xall




#belong to parent class: 
#getParameters01, diagonalizeMatrix, find_x_coord, setGraphLimits

#belong to child class, i.e. specific to BC and oscillationtype(transverse of longitudinal): 
#ploAnimation, setMatrix, sineCurves, also maybe setGraphLimits

#time defined in find_x_coord is repeated in plotAnimation - make a method for setting the time
#same as above goes for setting table length and equilibrium positions of masses.

#################################################################################################################
#################################################################################################################

class longFixed(Oscillator):
    def __init__(self,
                 numberMasses = 5,
                 relAmp = [1.0 , 0.8 , 0.6 , 0.4 , 0.4],
                 howLong = 4.0*np.pi*5.0,
                 howMany = 4000):  
                 #in the above YOU specify the arguments to create an instance of the "longFixed" class
        Oscillator.__init__(self,
                            numberMasses,
                            relAmp,
                            howLong,
                            howMany) 
                            # Now the parent class "Oscillator" is instanced with the parameters
                            # you specified in the previous line for the "longFixed" class - Hence
                            # both "Oscillator" and "longFixed" instances have the same parameters


    def getLengthTableLF(self):
        return float(self.numberMasses + 1)

    def getEqbPositionsLF(self):
        tabLength = self.getLengthTableLF()
        springLength = tabLength/float(self.numberMasses + 1)
        N = self.numberMasses
        x0 = []
        for ii in range(N):
            x0.append((-tabLength*0.5)+(ii+1)*springLength)
        return x0

    def setGraphLimitsLF(self):
        N = self.numberMasses
        xmin = -(float(N)+1)/2.0
        xmax = -xmin
        ymin = xmin
        ymax = xmax
        return (xmin,xmax,ymin,ymax)



    def sineCurvesLF(self):
        pi = np.pi
        N = self.numberMasses
        xmin,xmax,ymin,ymax = self.setGraphLimitsLF()
        sine_size = 50.0*float(N)
        x_sine_curve = np.arange(xmin,xmax,(xmax-xmin)/float(sine_size))
        y_sine_curve = []
        for ii in range(N):
            lmbda = 2.0*(xmax-xmin)/float(ii+1)
            kx = 2.0*pi/lmbda
            y_sine_curve.append(0.5*np.sin(kx*(x_sine_curve-xmin)))
        return x_sine_curve , y_sine_curve


    def plotAnimationLF(self):
        N = self.numberMasses
        Dn = self.setMatrix("fixed")
        (freq,amp) = self.diagonalizeMatrix(Dn,"fixed")
        xall = self.find_x_coord(freq,amp,self.getEqbPositionsLF())
        nmax = 0.5*(float(N)-1.0)
        nmin = -0.5*(float(N)-1.0)
        bias = np.arange(nmax,nmin-1,-1)
        time = self.setTime() #time = np.linspace(0 , self.howLong , self.howMany)
        timeSize = time.size

        #this was repeated in find_x_coord:
        tabLength = float(N)+1.0              # length of the table
        springLength = tabLength/(float(N+1)) # unstretched length of spring
        x0 = []
        for ii in range(N):
            x0.append((-tabLength*0.5)+(ii+1)*springLength) #equilibrium positions of the masses

        nspr = 20
        dy = 0.06

        xmin,xmax,ymin,ymax = self.setGraphLimitsLF()
        x_sine_curve , y_sine_curve = self.sineCurvesLF()

        fig = plt.figure()
        ax = plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
        ax.set_aspect("equal")
        ax.text(xmin,ymax+0.1,'$\copyright$2017, Bhaskar Kamble',fontsize=10,rotation=0)
        ax.axis("off")


        lines = []

        for ii in range(N):
            lobj = ax.plot([],[],lw=1,color="black")[0] # sine curves
            lines.append(lobj)

        for ii in range(N):
            lobj = ax.plot([],[],lw=1,color="black",linestyle = "-.")[0] # vertical dotted lines
            lines.append(lobj)

        for ii in range(N):
            for jj in range(N+1):
                lobj = ax.plot([],[],lw=1,color="black")[0]
                lines.append(lobj)
        for ii in range(N):
            for jj in range(N):
                lobj = ax.plot([],[],linestyle="none",marker="s",color="red",markersize=8)[0]
                lines.append(lobj)

        def init():
            for line in lines:
                line.set_data([],[])
            return lines


        def animate(i):
            #counter = 0
            for lnum,line in enumerate(lines):
                for ii in range(N):
                    if lnum == ii:
                        line.set_data(x_sine_curve,y_sine_curve[ii]+bias[ii])
                for ii in range(N):
                    if lnum == N+ii:
                        xeqb = np.array([x0[ii],x0[ii]])
                        yeqb = np.array([ymin,ymax])
                        line.set_data(xeqb,yeqb)
                for ii in range(N):
                    for jj in range(N+1):
                        if jj==0:
                            sp1 = xmin
                            sp2 = xall[N*ii][i]
                        if jj==N:
                            sp1 = xall[N*ii+N-1][i]
                            sp2 = xmax
                        if (jj != 0) and (jj != N):
                            sp1 = xall[N*ii+jj-1][i]
                            sp2 = xall[N*ii+jj][i]
                        if lnum == 2*N + ii*(N+1) + jj:
                            xspring = np.linspace(sp1,sp2,nspr)
                            yspring_temp = np.arange(xspring.size)
                            yspring = bias[ii] + dy*((-1.0)**yspring_temp)
                            line.set_data(xspring,yspring)
                for ii in range(N):
                    for jj in range(N):
                        if lnum == 2*N + N*(N+1)+N*ii+jj:
                            xls = xall[N*ii+jj][i]
                            yls = bias[ii]
                            line.set_data(xls,yls)
            return lines

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=timeSize, interval=50, blit=True)
        plt.show()













#################################################################################################################
#################################################################################################################

class longOpen(Oscillator):
    def __init__(self,
                 numberMasses = 5,
                 relAmp = [0.0 , 0.4 , 0.4 , 0.4 , 0.4],
                 howLong = 4.0*np.pi*5.0,
                 howMany = 4000):  
                 #in the above YOU specify the arguments to create an instance of the "longFixed" class
        Oscillator.__init__(self,
                            numberMasses,
                            relAmp,
                            howLong,
                            howMany) 
                            # Now the parent class "Oscillator" is instanced with the parameters
                            # you specified in the previous line for the "longFixed" class - Hence
                            # both "Oscillator" and "longFixed" instances have the same parameters



    def getLengthTableLO(self):
        return float(self.numberMasses)

    def getEqbPositionsLO(self):
        tabLength = self.getLengthTableLO()
        springLength = tabLength/float(self.numberMasses)
        N = self.numberMasses
        x0 = []
        for ii in range(N):
            #x0.append((-tabLength*0.5)+(ii+1)*springLength) #equilibrium positions of the masses
            x0.append((-tabLength*0.5 + 0.5*springLength)+ii*springLength)
        return x0


    def setGraphLimitsLO(self):
        N = self.numberMasses
        xmin = -float(N)/2.0
        xmax = -xmin
        ymin = xmin
        ymax = xmax
        return (xmin,xmax,ymin,ymax)

    def sineCurvesLO(self): # MODIFY THIS FOR OPEN BC # YES TRIED...
        pi = np.pi
        N = self.numberMasses
        xmin,xmax,ymin,ymax = self.setGraphLimitsLO()
        #xmin = -float(N)/2.0
        #xmax = -xmin
        #ymin = xmin
        #ymax = xmax
        sine_size = 50.0*float(N)
        x_sine_curve = np.arange(xmin,xmax,(xmax-xmin)/float(sine_size))
        y_sine_curve = []
        for ii in range(N-1):  #CHANGED CF. FIXEDBC
            lmbda = 2.0*(xmax-xmin)/float(ii+1)                    #CHANGED CF. FIXEDBC
            kx = 2.0*pi/lmbda
            y_sine_curve.append(0.5*np.cos(kx*(x_sine_curve-xmin)))#CHANGED CF. FIXEDBC
        return x_sine_curve , y_sine_curve


# CHECKPINT: IM ON THIS COMMAND IN THE BELOW PLOTANIMATION FUNCTION: (freq,amp) = self.diagonalizeMatrix(Dn,"open")
# ACHTUNG: table length, and equilibrium positions are also specific to the BC!!!
# xmin, xmax, ymin, ymax are also specific to the BC!!!

    def plotAnimationLO(self):
        N = self.numberMasses
        Dn = self.setMatrix("open")
        (freq,amp) = self.diagonalizeMatrix(Dn,"open") #YOU HAVE TO OMIT THE SMALLEST EVALUE AND EVECTOR
        #freq = freq[1:] #eliminates the 1st evalue
        #amp = amp[1:]   #eliminates the 1st array of amplitudes
        xall = self.find_x_coord(freq,amp,self.getEqbPositionsLO())
        xall = xall[N:]
        nmax = 0.5*(float(N-1)-1.0) #changed for open
        nmin = -0.5*(float(N-1)-1.0)#changed for open
        bias = np.arange(nmax,nmin-1,-1)
        time = self.setTime()
        timeSize = time.size

        #this was repeated in find_x_coord:
        tabLength = self.getLengthTableLO() #float(N)                  # length of the table
        springLength = tabLength/(float(N)) # unstretched length of spring
        x0 = self.getEqbPositionsLO()#[]
        #for ii in range(N):
        #    x0.append((-tabLength*0.5)+(ii+1)*springLength) #equilibrium positions of the masses

        nspr = 20
        dy = 0.06

        xmin,xmax,ymin,ymax = self.setGraphLimitsLO()
        #xmin = -float(N)/2.0
        #xmax = -xmin
        #ymin = xmin
        #ymax = xmax
        x_sine_curve , y_sine_curve = self.sineCurvesLO()

        fig = plt.figure()
        ax = plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
        ax.set_aspect("equal")
        ax.text(xmin,ymax+0.1,'$\copyright$2017, Bhaskar Kamble',fontsize=10,rotation=0)
        ax.axis("off")


        lines = []

        for ii in range(N-1):
            lobj = ax.plot([],[],lw=1,color="black")[0] # sine curves
            lines.append(lobj)

        for ii in range(N):
            lobj = ax.plot([],[],lw=1,color="black",linestyle = "-.")[0] # vertical dotted lines
            lines.append(lobj)

        for ii in range(N-1):
            for jj in range(N-1):
                lobj = ax.plot([],[],lw=1,color="black")[0] #springs
                lines.append(lobj)
        for ii in range(N-1):
            for jj in range(N):
                lobj = ax.plot([],[],linestyle="none",marker="s",color="red",markersize=8)[0] #masses
                lines.append(lobj)

        def init():
            for line in lines:
                line.set_data([],[])
            return lines


        def animate(i):
        #counter = 0
            for lnum,line in enumerate(lines):
                for ii in range(N-1): #sine curves
                    if lnum == ii:
                        line.set_data(x_sine_curve,y_sine_curve[ii]+bias[ii])
                for ii in range(N):   #vertical dotted lines
                    if lnum == N-1+ii:
                        xeqb = np.array([x0[ii],x0[ii]])
                        yeqb = np.array([ymin,ymax])
                        line.set_data(xeqb,yeqb)
                for ii in range(N-1):
                    for jj in range(N-1):
                        sp1 = xall[N*ii+jj][i]
                        sp2 = xall[N*ii+jj+1][i]
                        if lnum == (2*N-1) + ii*(N-1) + jj:
                            xspring = np.linspace(sp1,sp2,nspr)
                            yspring_temp = np.arange(xspring.size)
                            yspring = bias[ii] + dy*((-1.0)**yspring_temp)
                            line.set_data(xspring,yspring)
                for ii in range(N-1):
                    for jj in range(N):
                        if lnum == (2*N-1) + (N-1)*(N-1)+N*ii+jj:
                            xls = xall[N*ii+jj][i]
                            yls = bias[ii]
                            line.set_data(xls,yls)
            return lines

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=timeSize, interval=50, blit=True)
        plt.show()















#osc1 = longOpen()
#osc1.plotAnimation()



        #code
        #code
        


# the following classes are further possible:

#class longClosed(Oscillator):


#class transOpen(Oscillator):


#class transFixed(Oscillator):

