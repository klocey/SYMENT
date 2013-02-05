#! /usr/local/bin/python

import numpy as np
import math
import random

#FIND AND OPEN THE STARTING GENOME FILE
FILE1 = open('C:\Documents and Settings\icarus\My Documents\genome.txt','r') # Access a directory with a fasta file, the fasta file will be used to populate the community
genm = FILE1.read()
genm = genm.rstrip('\n')
s = len(genm)# Size of starting genome
#------------------------------------------
#THESE FOUR VARIABLES ARE THE ONLY ONES YOU NEED CHANGE
COMM = [400]   # Size of the community, the product of two equal sized dimensions
GENS = [1000]  # Number of generations
EXP = [5]      # The positive exponent of the mutation rate per base per individual reproductive event
edge = 0       # Number of rows and columns that will be removed to reduce edge effects
DISP = [0.5]   # Probability value for the geometic distribution; used to derive a dispersal kernal. 
#------------------------------------------

#MASTER LOOP: LOOPS THROUGH COMMUNITY SIZE (J), GENERATIONS (i.e. TIME), MUTATION RATES (EXP), and Dispersal kernels (DISP)
for J in COMM:
    for gens in GENS:
        for exp in EXP:
            for disp in DISP:
#------------------------------------------
# VARIABLES
                time = J*gens     # Number time steps to run
                comm = []         # An array to hold the community (an n by n square landscape of genomes)
                D = math.sqrt(J)  # Dimension of the community
                max1 = D - 1 - edge  # Max for rows and columns of the core
                min1 = edge        # Min for rows and columns of the core
                c = 0
                n = 0
                pm = 0
                x = 0
                row = 0
                colm = 0

#------------------------------------------
#POPULATE COMMUNITY WITH IDENTICAL GENOMES
                while n < J:# a loop to generate copies of genomes for spatially explicit community
                    g = list()
                    if colm == (D - 1): # Move to the next row when the last column is reached
                        g.append(row)  # Each element of the community array is an array containing a row coordinate
                        g.append(colm) # a column coordinate
                        g.append(genm) # and a genome
                        comm.append(g) # append the genome and its coordinates to the community array
                        row += 1    # move to next row
                        colm = 0    # reset column to 0
                        n += 1      # move to next element in community array
                    else:
                        g.append(row) 
                        g.append(colm)
                        g.append(genm)
                        comm.append(g)
                        colm += 1   # move to next column
                        n += 1      # move to next element in community array
                FILE1.close() # A community array of indentical genomes and their spatial coordinates is now created
#------------------------------------------
#INDUCE RANDOM REPLACEMENT OF BASES AND GENOMES OVER TIME
                while c < time: # Iterate until time is up.
                    ctdown = time - c
                    #print 'countdown = ',ctdown
                    x = random.randint(0,J-1)  # A uniform random number, used to randomly pick a genome from the community
                    genm = comm[x][2] # A randomly picked genome
                    row = comm[x][0]  # Row of the randomly picked genome
                    colm = comm[x][1] # Column of the randomly picked genome
                    rb = np.random.binomial(J*s,10**(-1*exp),1) # will a point mutation occur in this reproductive (i.e. replacement) event?    
                    if rb[0] == 1:    # PERFORM POINT MUTATION
                        b = random.randint(0,3)    # A uniform random number, used to randomly pick a replacing base (A,C,T, or G)  		
                        if b == 0:
                            pm = 'A'
                        elif b == 1:
                            pm = 'G'
                        elif b == 2:
                            pm = 'C'
                        elif b == 3:
                            pm = 'T'
                        r = random.randint(0,s-1)
                        g = list(genm)
                        g[r] = str(pm)
                        genm = "".join(g) # performing the point mutation on a random position within the genome

# PERFORM GENOME REPLACEMENT
                    rdist = int(np.random.geometric(disp,1)) # Distance moved across row
                    rdir = random.randint(1,2)  # Direction moved across row

                    cdist = int(np.random.geometric(disp,1)) # Distance moved across column
                    cdir = random.randint(1,2)  # Direction moved across column

                    if rdir == 1:   # Decrease row
                        if row - rdist >= min1: # if dispersal distance does not excede lower row bound
                            row -= rdist
                        else:
                            row = min1   # if disperal distance excedes lower row bound, move to the lower bound

                    elif rdir == 2: # Increase row
                        if row + rdist <= max1: # if dispersal distance does not excede upper row bound
                            row += rdist
                        else:
                            row = max1 # if disperal distance excedes upper row bound, move to the upper bound

                    if cdir == 1:   # Decrease column
                        if colm - cdist >= min1: # if dispersal distance does not excede lower column bound
                            colm -= cdist
                        else:
                            colm = min1 # if disperal distance excedes lower column bound, move to the lower bound

                    elif cdir == 2: # increase column
                        if colm + cdist <= max1: # if dispersal distance does not excede upper column bound
                            colm += cdist 
                        else:
                            colm = max1 # if dispersal distance excedes upper column bound, move to the upper bound

                    pos = int(row*D + colm)# For the genome to be replaced, derive its array index from the row & column coordinates
                    comm[pos][2] = genm  #...then replace its genome
                    c += 1
#PRINT THE COMMUNITY TO QUERRY AND REFERENCE FILES. THE COMMUNITY BECOMES A MULTI-FASTA FILE WHERE THE NAME OF THE GENOME IS THE ARRAY AND SPATIAL COORDINATES
                b = 0
                n = 0
                while n < J: # disregarding genomes on the edges
                    if comm[n][0] < min1 or comm[n][0] > max1 or comm[n][1] < min1 or comm[n][1] > max1: 
                        comm.pop(n)
                    n += 1

                M = len(comm)
                b = 0
                J = str(J)      #convert to string, to use in filehandle
                gens = str(gens)#convert to string, to use in filehandle
                exp = str(exp)  #convert to string, to use in filehandle

                disp = str(disp)
                qfile = 'C:\Documents and Settings\icarus\My Documents\Fasta\Q-'+J+'comm-'+gens+'gens-'+exp+'exp-'+disp+'disp.fasta'

                QFILE = open(qfile,'w') 

                while b < M:
                    row  = comm[b][0]
                    colm = comm[b][1]
                    genm = comm[b][2]
                    # Print the genome, its array coordinate, and row and column coordinates to the files created above
                    print >>QFILE, '>Q',b,row,colm,'\n',genm,'\n'

                    b += 1

                QFILE.close()
                J = int(J)
                gens = int(gens)
                exp = int(exp)
                #END REPLACEMENT LOOP
            #END DIST LOOP
        #END EXP LOOP  
    #END GENS LOOP
#END LOOP
print 'DONE'
