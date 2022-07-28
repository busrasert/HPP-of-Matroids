#############################################################################
##Authors: 
##Mario Kummer-Buesra Sert
##E-mail Addresses:
##mario.kummer@tu-dresden.de
##buesra.sert@tu-dresden.de
#############################################################################

##Catalog of Matroids from:
##https://doc.sagemath.org/html/en/reference/matroids/sage/matroids/matroids_catalog
Cat=[matroids.named_matroids.AG23minus(),matroids.named_matroids.AG32prime(),
    matroids.named_matroids.BetsyRoss(),matroids.named_matroids.Block_9_4(),
    matroids.named_matroids.Block_10_5(),matroids.named_matroids.D16(),
    matroids.named_matroids.ExtendedBinaryGolayCode(),
    matroids.named_matroids.ExtendedTernaryGolayCode(),
    matroids.named_matroids.F8(),
    matroids.named_matroids.Fano(),matroids.named_matroids.J(),
    matroids.named_matroids.K33dual(),matroids.named_matroids.L8(),
    matroids.named_matroids.N1(),matroids.named_matroids.N2(),
    matroids.named_matroids.NonFano(),matroids.named_matroids.NonPappus(),
    matroids.named_matroids.NonVamos(),matroids.named_matroids.NotP8(),
    matroids.named_matroids.O7(),matroids.named_matroids.P6(),
    matroids.named_matroids.P7(),
    matroids.named_matroids.P8(),matroids.named_matroids.P8pp(),
    matroids.named_matroids.P9(),
    matroids.named_matroids.Pappus(),matroids.named_matroids.Q6(),
    matroids.named_matroids.Q8(),
    matroids.named_matroids.Q10(),matroids.named_matroids.R6(),
    matroids.named_matroids.R8(),
    matroids.named_matroids.R9A(),matroids.named_matroids.R9B(),
    matroids.named_matroids.R10(),
    matroids.named_matroids.R12(),matroids.named_matroids.S8(),
    matroids.named_matroids.T8(),
    matroids.named_matroids.T12(),matroids.named_matroids.TernaryDowling3(),
    matroids.named_matroids.Terrahawk(),
    matroids.named_matroids.TicTacToe(),matroids.named_matroids.Vamos(),
    matroids.Wheel(4),matroids.Whirl(4),matroids.AG(2, 3)]

##We define functions to keep track of the name of a matroid from the catalog
##while adding more matroids to the list
    
def Dual(M):##Function for keeping the name of a matroid M while taking its dual
     DM=M.dual()
     RepM=repr(M)##Name of the matroid from the catalog
     r=RepM.rfind(":")
     DM.rename("Dual("+RepM[:r]+"): "+repr(DM)) ##Naming the dual
     return (DM)

CatD=[Dual(m) for m in Cat ]##List of duals of matroids from the catalog

    
def Ext(M): ##Function for keeping the name of a matroid M in its extension
     ExtM=M.extension() ##Extension of the matroid
     RepM=repr(M) ##Name of the matroid from the catalog
     r=RepM.rfind(":")
     ExtM.rename("Ext("+RepM[:r]+"): "+repr(ExtM)) ##Naming the extension
     return (ExtM)

ExtCat=[Ext(m) for m in Cat]##List of extension of matroids form the catalog 
   

def CoExt(M): ##Function for keeping the name of a matroid M in its co-extension
     CoExtM=M.coextension() ##Co-extension of the matroid
     RepM=repr(M) ##Name of the matroid from the catalog
     r=RepM.rfind(":")
     CoExtM.rename("CoExt("+RepM[:r]+"): "+repr(CoExtM)) ##Naming the co-extension
     return (CoExtM)

CoextCat=[CoExt(m) for m in Cat]##List of co-extension of matroids form the catalog

##Putting all the lists together
CAT=Cat+CatD+ExtCat+CoextCat ##Big catalog



% macaulay2--Runs Macaulay2 in Sage
----> Switching to Macaulay2 <--  
f="/home/.../.../BofHpcand"; --Calls the file: fill in the file location
--Warning: This is the list of 309 matroids, NOT the list of all matroids on 8 elements   
BofHpcand=value get f; --It calls the list in Macaulay2 syntax
--The syntax is of the form {{set{},...,set{}},...,{set{},...,set{}}}
--Ctrl-D
---> Exiting back to Sage <--

def convert(s): ## Function to convert M2 syntax to Sage syntax for list of Bases
    conv=list(eval(s.replace("{","[").replace("set","").replace("}","]")[1:-1])) 
    return (conv)
     
BofMs=convert(str(macaulay2("BofHpcand")))##List of Bases in Sage syntax
##The syntax is of the form [[[],...,[]],...,[[],...,[]]]


##Indices of 22 matroids that don't have Hpp in the list Hpcand
NotHpp=[73,76,77,81,82,83,85,88,92,93,94,95,96,97,99,100,101,112,113,114,116,117]

##Isomorphism Test on NotHpp
for i in NotHpp:
     B=BofMs[i]
     M=Matroid(bases=B)
     Iso=[]
     for j in range(len(CAT)):
          Mcat=CAT[j]
          if M.is_isomorphic(Mcat):
               Iso.append(Mcat)
     print("-------------------")
     print(i,Iso)     

##Testing which minors they have
for i in NotHpp:
     B=BofMs[i]
     M=Matroid(bases=B)
     Hasminor=[]
     for j in range(len(CAT)):
          Mcat=CAT[j]
          E=Mcat.groundset()
          n=len(E)
          if n<9:
               if M.has_minor(Mcat):
                    Hasminor.append(Mcat)
     print("----------------------")
     print(i,Hasminor)

##Testing which matroids have them as minors
##Loop takes a lot of time, we recommend to do them one by one
i=NotHpp[0]
B=BofMs[i]
M=Matroid(bases=B)
IsminorOf=[]
for j in range(len(CAT)):
     Mcat=CAT[j]
     E=Mcat.groundset()
     n=len(E)
     if n>8:
          if Mcat.has_minor(M):
               IsminorOf.append(Mcat)
print(i,IsminorOf)
     
     
##Indices of 14 matroids that don't pass SosDet test in the list Hpcand
NotDet=[58, 62, 63, 64, 66, 67, 68, 84, 86, 90, 106, 108, 109,110]

##Isomorphism Test on NotDet
for i in NotDet:
     B=BofMs[i]
     M=Matroid(bases=B)
     Iso=[]
     for j in range(len(CAT)):
          Mcat=CAT[j]
          if M.is_isomorphic(Mcat):
               Iso.append(Mcat)
     print("-------------------")
     print(i,Iso)
##Hpcand_110 is isomorphic to the Vamos matroid.
##None of the other matroids are isomorphic to a matroid from CAT

