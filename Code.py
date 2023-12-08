import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import csv
import copy

empty = 0
per = True
def PrintMat(array):
    str_ = ""
    for i in range(len(array)):
        for j in range(len(array)):
            str_ += str(array[i][j])
            str_ += " "
        print(str_)
        str_ = ""

def OneDimCase(size, self_e, bond):
    hamil = []
    
    row = []
    for i in range(size):
        for j in range(size):
            row.append(0)
        hamil.append(row)
        row = []
    
    for i in range(size):
        hamil[i][i] = self_e;
        
    
    for i in range(size - 1):
        hamil[i][i + 1] = bond;
        
        hamil[i+1][i] = bond;
        
    
    if(per == True):
        hamil[0][size - 1] = bond;
        hamil[size - 1][0] = bond;
    
    return hamil;

def TwoDimCase(size, self_e, bond):
    row = []
    col = []
    t_size = size*size
    for i in range(t_size):
        for j in range(t_size):
            if i == j:
                row.append(self_e)
            elif j == i +1 and (j + 1) % size != 1:
                row.append(bond)
            
            elif (j == i - 1 and (j + 1) % size != 0):
                row.append(bond)
            
            elif j % size == i % size and (j - size == i):
                row.append(bond)
            
            elif j % size == i % size and i - size == j:
                row.append(bond)
            
            elif((j + size) % t_size == i):
                row.append(bond)
            
            elif((i + size) % t_size == j):
                row.append(bond)
            
            else:
                row.append(empty)
        col.append(row)
        row = []
        
    for k in range(size):
        col[size*k][size*k + size - 1] = bond
        col[size*k + size - 1][size*k] = bond
    
        
        
        
        
        
    ham = np.array(col)
    return ham

def ThreeDimCase(size, self_e, bond):
    
    per = True
    
    len = size*size*size;
    sheet_size = size*size;
    hamil = []
    row = []
    
    #for i in range(len):
    for j in range(len):
        row.append(0)
    for i in range(len):
        hamil.append(copy.deepcopy(row))
        
    
    i = 0
    for i in range(len):
        line = hamil[i]
        line[i] = self_e
        
    
    
    for i in range(len - 1):
        if((i+1) % size != 0):
            hamil[i][i + 1] = bond;
            hamil[i+1][i] = bond; #//conjugate here
        
    
    for i in range(len - sheet_size):
        hamil[i][sheet_size + i] = bond;
        hamil[sheet_size + i][i] = bond; #// conjugate
    

    for i in range(len - size):
        if(i % sheet_size < sheet_size - size):
            hamil[i + size][i] = bond;
            hamil[i][ i + size] = bond; #// conjugate 
        
    
    
    
    if(per == True):

        for i in range(len - 2*size):
            if(i % sheet_size < size):
                hamil[i][i + sheet_size - size] = bond;
                hamil[i + sheet_size - size][i] = bond; #// conjugate
            
        

        for i in range(sheet_size):
            hamil[i * size][(i+1)*size - 1] = bond;
            hamil[(i+1)*size - 1][i * size] = bond; #//conjugate
            
            hamil[i][i + len - sheet_size] = bond;
            hamil[len - sheet_size + i][i] = bond; #// conjugate
        
    return hamil

def OneDimSpin(size, self_e, bond):
    
    hamil = []
    row = []
    self_up = self_e
    self_down = -1*self_e
    
    len = size*2
    
    for i in range(len):
        for j in range(len):
            row.append(0)
        hamil.append(row)
        row = []
    
    
    i = 0
    for i in range(size):
        hamil[i][i] = self_up;
        hamil[i+size][i+size] = self_down;
    
    for i in range(size - 1):
        hamil[i][ i + 1] = bond;
        hamil[i + size][i + 1 + size] = bond;
        hamil[i+1][i] = bond;                  #//Add conjugate
        hamil[i+1 + size][i + size] =bond;
    
    if per == True:
        hamil[0][size - 1] = bond;
        hamil[size - 1][0] = bond;
        hamil[size][size*2 - 1] = bond;
        hamil[size*2 - 1][size] = bond;
        
    return hamil

def TwoDimSpin(size, self_e, bond):
    self_up = self_e
    self_down = -1*self_e
    
    len = size*size;
    hamil = []
    row = []
    for i in range(2*len):
        for j in range(2*len):
            row.append(0)
        hamil.append(row)
        row = []

    for i in range(len):
       hamil[i][ i] = self_up;
       hamil[i + len][i + len] = self_down;
        
    
    for i in range(len - 1):
        if((i+1) % size != 0):
            hamil[i][i + 1] = bond;
            hamil[i+1][ i] = bond; #//conjugate here

            hamil[i + len][i + 1 + len] = bond;
            hamil[i + 1 + len][ i + len] = bond; #// conjugate here
        
    
    for i in range(len - size):
        hamil[i + size][ i] = bond;
        hamil[i][i + size] = bond; #// conjugate here

        hamil[i + size + len][i + len] = bond;
        hamil[i + len][i + size + len] = bond;
    
    if per == True:
        for i in range(size):
            hamil[i][i + len - size] = bond;
            hamil[i + len - size][ i] = bond;# // conjugate

            hamil[i * size][(i+1)*size - 1] = bond;
            hamil[(i+1)*size - 1][i * size] = bond; #// conjugate

            hamil[i + len][ i + len - size + len] = bond;
            hamil[i + len - size + len][i + len] = bond; #// conjugate

            hamil[i * size + len][(i+1)*size - 1 + len] = bond;
            hamil[(i+1)*size - 1 + len][i * size + len] = bond; #//conjugate

        
    return hamil

def ThreeDimSpin(size, self_e, bond):
    self_up = self_e
    self_down = -1*self_e
    
    hamil = []
    row = []
    len = size*size*size;
    sheet_size = size*size;
    
    for i in range(2*len):
        for j in range(2*len):
            row.append(0)
        hamil.append(row)
        row = []

    for i in range(len):
        hamil[i][i] = self_up;
        hamil[i+len][ i+len] = self_down;
        
    
    for i in range(len - 1):
        if((i+1) % size != 0):
            hamil[i][ i + 1] = bond;
            hamil[i+1][i] = bond;

            hamil[i + len][i + 1 + len] = bond;
            hamil[i + 1 + len][i + len] = bond;
        
    
    for i in range(len - sheet_size):
        hamil[i][ sheet_size + i] = bond;
        hamil[sheet_size + i][ i] = bond;
        
        hamil[i + len][sheet_size + i + len] = bond;
        hamil[sheet_size + i + len][i + len] = bond;
    

    for i in range(len - size):
        if(i % sheet_size < sheet_size - size):
            hamil[i + size][ i] = bond;
            hamil[i][i + size] = bond;

            hamil[i + size + len][i + len] = bond;
            hamil[i + len][i + size + len] = bond;
        
    
    
    
    if(per ==True):
        for i in range(len - 2*size):
            if(i % sheet_size < size):
                hamil[i][i + sheet_size - size] = bond;
                hamil[i + sheet_size - size][ i] = bond;

                hamil[i + len][ i + sheet_size - size + len] = bond;
                hamil[i + sheet_size - size + len][ i + len] = bond;
            
        

        for i in range(sheet_size):
            hamil[i * size][ (i+1)*size - 1] = bond;
            hamil[(i+1)*size - 1][ i * size] = bond;

            hamil[i * size + len][ (i+1)*size - 1 + len] = bond;
            hamil[(i+1)*size - 1 + len][ i * size + len] = bond;
            
            hamil[i][ i + len - sheet_size] = bond;
            hamil[len - sheet_size + i][ i] = bond;

            hamil[i + len][ i + len - sheet_size + len] = bond;
            hamil[len - sheet_size + len + i][ i + len] = bond;
    
    
    
    return hamil
    

def DenStat(E, eta): #E is a list of eigenvalues // assumes all real values

    mn=min(E)-0.1
    mx=max(E)+.1 #E eigenvalues
    Eo = []
    
    g1 = []
    i = mn
    
    while i < mx:
        Eo.append(i)
        i += 0.005
            
    i = 1
    ee = []
    j = 0

    while j < len(Eo):
        ee.append(j)
        
        GG1=[]
        for i in range(len(E)):
            x = E[i]
            #GG.append(1/3.14*eta/((Eo[j]-x)**2 +eta**2)/(1)) # lorentzian
    
            GG1.append(1/(eta*(2*3.14)**(1/2))*2.72**(-(Eo[j]-x)**2/(2*eta**2))/(1)) # gaussian
        
        
        g1.append(np.sum(GG1))
        j += 1
    
    
    #plt.plot(Eo, g, 'r--', Eo, g1, 'b')
    plt.plot(Eo, g1, 'b')
    plt.show()
    
    return g1, Eo
    
def IntegStates(f_x):

    delta = 0.005;
    integral = 0;
    step = 0;

    coef = (delta/2);

    integral += f_x[0]*coef;
    integral += f_x[-1]*coef;

    for i in range(len(f_x) - 1):
        step += f_x[i];
    
    integral += step*delta;
    return integral;

    
def CSV_Writer(ham):
    fp = open("output.dat", "w")
    writer = csv.writer(fp)
    
    for i in range(len(ham)):
        writer.writerow(ham[i])
    fp.close()
        
def main():
    size = int(input('Enter the number of lattice sites in a row: '))
    dim = input('(Case sensitive) Is your lattice 1D, 2D, or 3D? ')
    spin = input('Would you like to model spin? (y/n): ')
    self_e = float(input('What is the on-site energy? '))
    bond = float(input('Enter the hopping parameter: '))
    per = input('Would you like to model periodic boundary conditions? (y/n): ')
    eta = float(input('What is your smoothing parameter? '))
    
    if per.lower() == 'n':
        per = False
    
    hamil = []
    
    
    if spin.lower() == 'n':
        if dim == '1D':
            hamil = OneDimCase(size, self_e, bond)
        elif dim == '2D':
            hamil = TwoDimCase(size, self_e, bond)
        elif dim == '3D':
            hamil = ThreeDimCase(size, self_e, bond)
        else:
            print('Oops! Something went wrong with these inputs.')
    elif spin.lower() == 'y':
        if dim == '1D':
            hamil = OneDimSpin(size, self_e, bond)
        elif dim == '2D':
            hamil = TwoDimSpin(size, self_e, bond)
        elif dim == '3D':
            hamil = ThreeDimSpin(size, self_e, bond)
        else:
            print('Oops! Something went wrong with these inputs.')
    else:
        print('Oops!, Something went wrong with these inputs')
    hamil = np.array(hamil)
    eigs, eig = la.eigh(hamil)
    f_x, x = DenStat(eigs, eta)
    cont = input('Does this graph look good? (y/n): ')
    while cont.lower() != 'y':
        eta = float(input('Enter a new smoothing parameter: '))
        f_x, x = DenStat(eigs, eta)
        cont = input('Does this graph look good? (y/n): ')
    
    print('The total number of states is: ', IntegStates(f_x))
    
main()
