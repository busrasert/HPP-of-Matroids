---------------------------------------------------------------------------------
--Authors: 
--Mario Kummer-Buesra Sert
--E-mail Addresses:
--mario.kummer@tu-dresden.de
--buesra.sert@tu-dresden.de
---------------------------------------------------------------------------------


needsPackage("SumsOfSquares");
needsPackage("Matroids");
needsPackage("NumericalAlgebraicGeometry");
needsPackage("SemidefiniteProgramming");
  
------------------------------------------------------------------------------
--<<<<<<<<<<<<<<<<<<Tests for the HPP>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
------------------------------------------------------------------------------

L8=allMatroids 8;--List of all non-isomorphic matroids with 8 elements
-- #L8=1724
L={};--List of matroids on 8 elements with rank 3 or 4
for M in L8 do (
    if rank M==3 or rank M==4
       then L=append(L,M);
       );
--#L=1265
S={}; --List of simple matroids
for M in L do (
    if isSimple M
       then S=append(S,M);
       );
--#S=685
C={}; --List of simple and connected matroids
for M in S do (
    if #components M==1
       then C=append(C,M);
       );
--#C=659  


F=specificMatroid "fano"; -- $F$
NF=specificMatroid "nonfano"; -- $F^{-}$
F2=relaxation(NF,set{0,1,6}); --$F^{--}$
K4e=relaxation(F2,set{2,1,5});--$M(K4)+e$
preF3=relaxation(NF,set{0,2,4});
F3=relaxation(preF3,set{0,3,5});--$F^{-3}$

DF=dual F; --dual of $F$
DNF=dual NF; --dual of $F^{-}$
DF2=dual F2; --dual of $F^{--}$
DK4e=dual K4e; --dual of $M(K4)+e$
DF3=dual F3; --dual of $F^{-3}$


H={}; --List of matroids that don't have the forbidden minors of rank 3
for M in C do (
    if not hasMinor(M,F) and not hasMinor(M,NF) and not hasMinor(M,F2)
      and not hasMinor(M,F3) and not hasMinor(M,K4e)
       then H=append(H, M);
       );
      
 Hpcand={};-- List of matroids that don't have any of the forbidden minors
 --We recommend to save this list in a file in order to call back when needed again
for M in H do (
    if rank M==4 and not hasMinor(M,DF) and not hasMinor(M,DNF)
    and  not hasMinor(M,DF2) and not hasMinor(M,DF3) and not hasMinor(M,DK4e)
       then Hpcand=append(Hpcand, M);
    if rank M==3
       then Hpcand=append(Hpcand, M);
       );
--#Hpcand=309


R=QQ[x_0 .. x_7]; --Ring of the basis generating polynomials
      
BgP=method();--Function for finding the basis generating polynomial of a matroid
BgP(Matroid,Ring):= RingElement => (M,R) ->(
    h:=0_R; B:=bases(M);
    for i from 0 to (#B-1) do (
       L:=toList(B_i); N:=1;
       for j from 0 to (#L-1) do (
           N=N*x_(L_j); --monomials
           );
     h=h+N
        );
     return h
     );
 --Creating list of indices for Rayleigh difference 
S1=toList(0..6);
S2=toList(0..7);
A=set(S1)**set(S2);
I=toList A;
J={};--the list of distinct pairs of indices i,j
for i from 0 to (#I-1) do (
    a:=I_i;
    if a_0!=a_1
      then J=append(J,a)
      );

RDSos=method();-- Function for SoS test on Rayleigh differences for a matroid 
RDSos(Matroid,List,Ring):= Boolean =>(M,J,R)-> ( 
     h:=BgP(M,R); --Basis generating polynomial
     SS={}; 
     for a in J do(
         i:=a_0; j:=a_1;
         rayl=diff(x_i,h)*diff(x_j,h)-h*diff(x_i,diff(x_j,h));--Rayleigh difference
	 sol=solveSOS rayl; --Sos test
	 if status(sol)=="SDP solved, primal-dual feasible"
            then SS=append(SS,{i,j});
	    );
      if #SS==0
         then return false
      else
         return true
       );


unkwn={};--List of matroids that don't pass SoS test
for i from 0 to #Hpcand-1 do (
    print i; --to keep track of the process
    M=Hpcand_i;--Matroid
    if not RDSos(M,J,R)
       then unkwn=append(unkwn,M)
       );
--It takes time and memory.
--We recommend to run this code on several parts of Hpcand instead of all in once
--For example, do it for i from 0 to 199 first, then for i from 200 to 308
--#unkwn=22


--List of indices of 22 matroids from unkwn in Hpcand  
NotHpp={73,76,77,81,82,83,85,88,92,93,94,95,96,97,99,100,101,112,113,114,116,117}

Hpp={};--List of indices of 287 Matroids from Hpcand that have HPP
for i from 0 to #Hpcand-1 do (
    if not member(i,NotHpp)
       then Hpp=append(Hpp,i);
       );

RDGram=method();--Function that returns the list of rings of Gram matrices
RDGram(Matroid,List,Ring):= List =>(M,J,R) -> ( 
    h:=BgP(M,R); --Basis generating polynomial
    SS:={}; 
    for a in J do(
	i:=a_0; j:=a_1;
    	rayl=diff(x_i,h)*diff(x_j,h)-h*diff(x_i,diff(x_j,h));
	sol=solveSOS rayl; --Sos test
    	if status(sol)=="SDP solved, primal-dual feasible" 
	   then (G=sol#GramMatrix; SS=append(SS, ring G)) --Ring of the Gram matrix 
        else
	   continue
	     );
    return SS
    );

--Recording the data in a file
f="/home/.../.../Gram"; --Fill in the file location   
for i from 0 to #Hpp-1  do (
    print i; j=Hpp_i;
    M=Hpcand_j;
    K=RDGram(M);
    f<< toString({i, K})<< ","<< flush --Writing on the file
    );
f<<close;--Closing the file
    
--It takes time and memory.
--We recommend to run this code on several parts of Hpp instead of all in once


--{42,44,46,48,49,50,52,53,54,59,60,61,70,72,79,120,141} 
--List of indices of matroids in Hpcand for which some Gram matrices have floating point entries


-------------------------------------------------------------------------------------
--<<<<<<<<<<<<<<<<<<<<<<<<<<<Tests on Hyperbolicity>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
-------------------------------------------------------------------------------------

--Creating the list of vectors with 0-1 entries
K=toList(0..2^8-1);
D={};--List of all directions with 0-1 entries
--Its elements are {0,0,1,0,1,1,1,0} etc.
for i from 0 to 2^8-1 do (
    a=();
    for j from 0 to 7 do(
         if K_i & (2^j) == 0
             then a=append(a,0)
         else a=append(a,1)
         );
     D=append(D,toList(a));
     );
U=CC[t];--Ring of univariate polynomial $h_M(et-v)$
--(over complex numbers as we will check complex solutions)
univP=method();--Function to create the univariate polynomial $h_{M}(et-v)$
univP(Matroid,List,List,Ring):= RingElement => (M,e,v,U) -> (
    B:=bases(M); p:=0_U;
    for i from 0 to #B-1 do(
        s:=toList(B_i); c:=1;
        for j from 0 to #s-1 do (
            k:=s_j; m=(e_k*t-v_k); c=c*m;
	    ); 
	 p=p+c;
	 );
    return p
    );
  


 --For each matroid, creating a file with polynomials with non-real roots
    
for i from 0 to #unkwn-1 do(
    S={}; --List of tuples ((e,v),univP(M,e,v,U))
    M=unkwn_i;--Matroid from the list of 22 matroids
    for j from 0 to 2^8-1 do(
        for k from 0 to 2^8-1 do(
           h=univP(M,D_k,D_j,U); --Univariate polynomial
           if h!=0
              then (d=degree h; 
              if #d!=0
                 then (F={h}; sys:=solveSystem F;
                 n:=#sys-#(realPoints sys); --Checks the number of non-real roots
                 if n!=0
                    then S=append(S,((D_k,D_j),h))););
               );
            );
      y="/home/.../.../NonRealM"|i; --Creates a file (fill in the file location)  
      y<< toString(S)<< close; --Writes on the file
      );  



PD={};--List of random positive directions
for i from 0 to 2^8-1 do (
    b=();
    for j from 0 to 7 do (
        b=append(b,random(0.,100.))
         );
    PD=append(PD,toList(b));
    );

ND={};--List of random directions
for i from 0 to 2^8-1 do (
    b=();
    for j from 0 to 7 do (
        b=append(b,random(-100.,100.))
         );
    ND=append(ND,toList(b));
    );

T={0,1,2,12,13,19}--List of indices of 7 matroids in unkwn
for i in T do(
    S={}; --List of tuples ((e,v),univP(M,e,v,U))
    M=unkwn_i;
    for j from 0 to 2^8-1 do(
        for k from 0 to 2^8-1 do(
           h=univP(M,PD_k,ND_j,U); --Univariate polynomial with random directions
           if h!=0
              then (d=degree h; 
              if #d!=0
                 then (F={h}; sys:=solveSystem F;
                 n:=#sys-#(realPoints sys); --Checks the number of non-real roots
                 if n!=0
                    then S=append(S,((PD_k,ND_j),h))););
               );
            );
      print S;--One can also print it in a file as above
      );  


-------------------------------------------------------------------------------------
--<<<<<<<<<<<<<<<<<<<<<Tests for Determinantal Representability>>>>>>>>>>>>>>>>>>>>
-------------------------------------------------------------------------------------


R=QQ[x_0 .. x_7]; --Ring of the basis generating polynomials
  
 --Creating list of indices for Rayleigh difference 
S1=toList(0..6);
S2=toList(0..7);
A=set(S1)**set(S2);
I=toList A;
J={};--the list of distinct pairs of indices i,j
for i from 0 to (#I-1) do (
    a:=I_i;
    if a_0!=a_1
      then J=append(J,a)
      );
      
SosDet=method();--Function for testing determinantal representability
SosDet(Matroid,List, Ring):= Boolean =>(M,J,R) -> ( 
     h:=BgP(M,R);--Basis generating polynomial 
     S:={}; 
     for a in J do(
         i:=a_0; j:=a_1;
	 rayl=diff(x_i,h)*diff(x_j,h)-h*diff(x_i,diff(x_j,h));--Rayleigh difference
	 sol=solveSOS rayl; --Sos test
	 if status(sol)=="SDP solved, primal-dual feasible"
            then S=append(S,{i,j});
	     );
      if #S!=#J and #S>0
         then return false--Implies that M doesn't have a determinantal representation
      else
         return true--Doesn't imply anything in terms of determinantal representation
         );
         
--Hpp is the list of indices of matroids with Hpp (defined above)
 
NotDet={};--List of indices of matroids that don't past SosDet test in Hpcand
--Test is applied on matroids with HPP 
for j from 0 to #Hpp-1 do (
    print j;--to keep track of the process
    i=Hpp_j;--Corresponding index in Hpcand
    M=Hpcand_i;  
    if not SosDet(M,J,R)
       then  NotDet=append(NotDet,i);
       );
--It takes time. We recommend to split indices into parts and run the test on the parts 
--NotDet={58, 62, 63, 64, 66, 67, 68, 84, 86, 90, 106, 108, 109,110}
--List of indices of matroids in Hpcand that don't have a determinantal
--representation

--Hpcand_110 is the Vamos Matroid

--L8IndexNotDet={409, 413, 414, 415, 417, 418, 419, 438, 440, 445, 498, 500, 501,502}
--List of indices of Matroids in L8 that don't have a determinantal representation


----------------------------------------------------------------------------------
--<<<<<<<<<<<<<<<<<<<<<<<<<<<Non-SOS Certificates>>>>>>>>>>>>>>>>>>>>>>>>>>>>
-----------------------------------------------------------------------------------

notSosRayl={409, 413, 414, 415, 417, 418, 419, 438, 440, 445, 498,
500, 501,502};--502 is Vamos
--List of L8 indices of Non-SOS Rayleigh matroids 

RDSosIndx=method();
-- Function for SoS test on Rayleigh differences and return
--the indices that make sos
RDSosIndx(Matroid,List,Ring):= List =>(M,J,R)-> ( 
     h:=BgP(M,R); --Basis generating polynomial
     SS={}; --List of indices
     for p in J do(
         i:=p_0; j:=p_1;
         rayl=diff(x_i,h)*diff(x_j,h)-h*diff(x_i,diff(x_j,h));
         --Rayleigh difference
	 sol=solveSOS rayl; --Sos test
	 if status(sol)=="SDP solved, primal-dual feasible"
            then SS=append(SS,(i,j));
	    );
      if #SS>0
         then return SS
       );

BgP=method();
--Function for finding the basis generating polynomial
BgP(Matroid,Ring):= RingElement => (M,R) ->(
    h:=0_R; B:=bases(M);
    for i from 0 to (#B-1) do (
       L:=toList(B_i); N=1;
       for j from 0 to (#L-1) do (
           N=N*x_(L_j); --monomials
           );
     h=h+N
        );
     return h
     );

L8=allMatroids 8;
M=L8_417;--A matroid from the list

F1=QQ[x_0..x_7];
--Polynomial ring of the basis generating polynomial

h=BgP(M,F1);--The basis generating polynomial

S1=toList(0..6);
S2=toList(0..7);
C=set(S1)**set(S2);
I=toList C;
J={};--the list of distinct pairs of indices i,j
for i from 0 to (#I-1) do (
    p:=I_i;
    if p_0!=p_1
      then J=append(J,p)
      );
  
N=RDSosIndx(M,J,F1);
--List of indices that make the Rayleigh difference Sos 
K={};
--List of the indices that fail to make Rayleigh difference Sos
for i in J do (
    if not member(i,N)
       then K=append(K,i);
);

i=K_0_0; j=K_0_1;--Setting indices i,j
rayl=diff(x_i,h)*diff(x_j,h)-h*diff(x_i,diff(x_j,h));
--The Rayleigh difference given by the first index set from K

n=binomial(6,3);--The number of monomials in the matrix X
k=n*(n+1)//2;--The number of distinct entries of G
RS=QQ[a_0..a_(k-1)][x_0..x_7];--Creating the ring
R=coefficientRing RS;
A=genericSymmetricMatrix(R,a_0,n); 
G=map(RS^n,RS^n,(i,j)-> (A_(i,j))_RS);
--For changing the ring of A
--This is the initial Gram Matrix with entries ai

y=toList(0..7)-set{i,j};
--List of indices of variables in the Rayleigh difference

substs=subsets(y,3);--Subsets of indices that can appear

E={};--List of monomials to insert in the matrix X
for i from 0 to #substs-1 do (
    c:=1;  s:=substs_i; 
    for j from 0 to #s-1 do(
        t:=s_j;  m=x_t;
        c=c*m; --Monomial
	);
    E=append(E,c)
    );
X=matrix{E};
--matrix M with entries that are prod. of x_i for each subset
XT=transpose X; 
Q=X*G*XT;--Matrix with one polynomial entry 

f=substitute(rayl,RS)-Q_(0,0);--Polynomial rayl-Q_00
Cf=coefficients f;
I=ideal(Cf);--Ideal of the relation of the entries
G=G%I;--Inserting relations in G
varS={};
for i from 0 to n-1 do (
    for j from 0 to n-1 do (
	gij=G_(i,j);--Entry Gij
	T=terms gij;
	for k in T do (
	    Lmk=leadMonomial substitute(k,R);
            --Lead monomial of the term
	    if member(Lmk,gens R)
	    	then
	    	  varS=append(varS,Lmk);
	    );
	);
    );

varS=toList(set(varS));--List of used variables (lambdas)

GiS={};--List of matrices Gi

--Creating matrices Gi
IR=ideal(gens R);
G0=substitute(G,R)%IR;--Matrix G_0
GiS=append(GiS,G0);

for i from 0 to #varS-1 do ( 
    C:=G-G0; 
    m:=varS_i;
    C=substitute(C,{m=>1});
    --Substitute 1 for the variable
    C=substitute(C,R)%IR;
    --Substitute 0 for the others
    GiS=append(GiS,C);
    );

Traces={};--List of traces of MGi

for i from 0 to #GiS-1 do (
    Gg=GiS_i; Z=A*Gg;
    if not member(trace Z, Traces) 
      then
       Traces=append(Traces,trace Z);
    );

Matrel={};--Matrices for the trace relation 

for i in Traces do(
    Mt=A;--Symmetric matrix with entries ai
    Ti=terms substitute(i,R);
    for j in Ti do(
	Lm=leadMonomial j;
	Lc=leadCoefficient j;
	Mt=substitute(Mt,{Lm=>Lc});
	--Inserting the coefficients from the trace
	);
    Mt=Mt%IR;--Setting other entries zero
    Matrel=append(Matrel,Mt);
    );


r=#GiS;
C=map(QQ^n,QQ^n,(i,j)->0);--We don't want to minimize something
b=map(QQ^r,QQ^1,(i,j)->0);--We want the trace of M*Gi=0

P=sdp(C,toSequence(GiS),b);--SDP
(X,y,Z,stat) = optimize P;
--X is a point from the solution set


Mvecs={};--List of vectorized matrices

for i from 0 to #Matrel-1 do(
    m=substitute(Matrel_i,QQ); Sv=smat2vec m;
    Mvecs=append(Mvecs,Sv);
    );


TMvecs=transpose matrix{Mvecs};
u=#Mvecs;
bb=map(QQ^u,QQ^1,(i,j)->0);

for i in GiS do (print trace(Mat*i))
(ispsd,Mat) = roundPSDmatrix(X,TMvecs,bb,10^5);
--Rounds the SDP solution, satisfying the trace condition
--ispsd
--Mat--Certification matrix

--for i in GiS do (print trace(Mat*i))
--Checking that all traces are zero

--ispsd
--Mat


-----------------------------------------------------------------------------------
--<<<<<<<<<<<<<<<<<<<<<<<<<<<<Non-Rayleigh Matroids>>>>>>>>>>>>>>>>>>>>>>>>>
-----------------------------------------------------------------------------------

M=L8_891;--Matroid from the list

R=QQ[x_0..x_7]
h=BgP(M,R);--Basis generating Polynomial

S1=toList(0..6);
S2=toList(0..7);
C=set(S1)**set(S2);
I=toList C;
J={};--the list of distinct pairs of indices i,j
for i from 0 to (#I-1) do (
    p:=I_i;
    if p_0!=p_1
      then J=append(J,p)
      );

K=J_27;--indices (i,j)
--For different L, insert the Jindx value given with L

I=K; k=I_0; t=I_1;
rayl=diff(x_k,h)*diff(x_t,h)-h*diff(x_k,diff(x_t,h));
--Rayleigh difference

Y=toList(set(toList(0..7))-set({k,t}));
--List of indices i of xi that appear in Rayl 
e0=Y_0; e1=Y_1; e2=Y_2; e3=Y_3; e4=Y_4; e5=Y_5;

L={1,19/100,141/100,141/100,19/100,19/100}--L8_891-Jindx27
--Point that gives the negatie value
--See below for list of L's

mini1=sub(rayl,{x_(e0)=>L_0, x_(e1)=>L_1, x_(e2)=>L_2,
 x_(e3)=>L_3, x_(e4)=>L_4, x_(e5)=>L_5});
print mini1;--prints the negative value


p={};--List of terms of h without xi and xj
q={};--List of terms of h with xi and xj
r={};--List of terms of h with xi without xj
s={};--List of terms of h With xj without xi

Trms=terms h;
Ex=exponents h;
for j from 0 to #Ex-1 do (
   I=toList(Ex_j);
   if I_k==0 and I_t==0
      --Without xi and xj
      then  
       p=append(p,Trms_j)
    );
for j from 0 to #Ex-1 do (
    I=Ex_j;
    if I_k==1 and I_t==1
       --With xi and xj
       then  
        q=append(q,Trms_j)
    );
for j from 0 to #Ex-1 do (
    I=Ex_j;
    if I_k==1 and I_t==0
       --With xi without xj
       then  
        r=append(r,Trms_j)
    );
for j from 0 to #Ex-1 do (
    I=Ex_j;
    if I_k==0 and I_t==1
       --With xj without xi
       then  
        s=append(s,Trms_j)
    );

polD=sum p;--polynomial d
polA=sum q;--polynomial a*xi*xj
polC=sum s;--polynomial c*xj
polB=sum r;--polynomial b*xi

--Substituting point L in polynomials a,b,c,d 
minid=sub(polD,{x_(e0)=>L_0, x_(e1)=>L_1, x_(e2)=>L_2,
 x_(e3)=>L_3, x_(e4)=>L_4, x_(e5)=>L_5});
--print minid;

minia=sub(polA,{x_(e0)=>L_0, x_(e1)=>L_1, x_(e2)=>L_2,
 x_(e3)=>L_3, x_(e4)=>L_4, x_(e5)=>L_5});
--print leadCoefficient minia;

minic=sub(polC,{x_(e0)=>L_0, x_(e1)=>L_1, x_(e2)=>L_2,
 x_(e3)=>L_3, x_(e4)=>L_4, x_(e5)=>L_5});
--print leadCoefficient minib;

minib=sub(polB,{x_(e0)=>L_0, x_(e1)=>L_1, x_(e2)=>L_2,
 x_(e3)=>L_3, x_(e4)=>L_4, x_(e5)=>L_5});
--print leadCoefficient minic;


dL=leadCoefficient minid;--d(L)
aL=leadCoefficient minia;--a(L)
bL=leadCoefficient minib;--b(L)
cL=leadCoefficient minic;--c(L)
--Checking if ad/bc>8/7
(dL*aL)/(bL*cL)>8/7


-------------------------------------------
--List of points for matroids on the table
-------------------------------------------

--L={1,2302/100000,2302/100000,88/100000,1936/100000,
--1936/100000}
--L8_768-Jindx27
---------------------------------

--L={1,4096/10000,4234/10000,493/10000,1071/10000,373/10000
--L8_816-Jindx2
---------------------------------

--L={1,9348/10000,4050/10000,4050/10000,799/10000,799/10000}
--L8_821-Jindx2
---------------------------------

--L={1,092/1000,1/10,7/1000,52/1000,4/1000}
--L8_825-Jindx2
---------------------------------

--L={1,2673/10000,2822/10000,735/10000,1128/10000,63/10000}
--L8_878-Jindx19
---------------------------------

--L={1,296/100,305/100,640/100,48/100,7/100}
--L8_879-Jindx19
---------------------------------

--L={1,94/1000,102/1000,11/1000,56/1000,2/1000}
--L8_882-Jindx19
---------------------------------

--L={1,13/100,13/100,2/100,4/100,4/1000}
--L8_883-Jindx19
---------------------------------

--L={1,19/100,141/100,141/100,19/100,19/100}
--L8_891-Jindx27
---------------------------------

--L={1,74/1000,73/1000,7/1000,73/1000,3/1000}
--L8_894-Jindx27
---------------------------------

--L={1,37/1000,36/1000,2/1000,36/1000,1/1000}
--L8_895-Jindx27
---------------------------------

--L={1,3/10000,3/10000,3/10000,3/10000,3/10000}
--L8_896-Jindx27
---------------------------------

--L={1,3/100,3/100,2/1000,3/100,1/1000}
--L8_910-Jindx27
---------------------------------

--L={1,2/100,2/100,1/1000,2/100,1/1000}
--L8_911-Jindx27
---------------------------------

--L={1,25/100000,25/100000,25/100000,25/100000,25/100000}
--L8_912-Jindx27
---------------------------------
