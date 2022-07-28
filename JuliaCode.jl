#######################################################
##Authors: 
##Mario Kummer-Buesra Sert
##E-mail Addresses:
##mario.kummer@tu-dresden.de
##buesra.sert@tu-dresden.de
######################################################

##Julia Code for Finding the Negative Points##


import Pkg
Pkg.add("HomotopyContinuation")
Pkg.add("IterTools")

using IterTools
subsets(s::AbstractSet, k) = IterTools.subsets(collect(s), k)

function revLex(f1, f2)
    if maximum(f1, dims=1) < maximum(f2,dims=1)
        return true;
    elseif maximum(f1,dims=1) > maximum(f2,dims=1) 
        return false;    
    else  
        if length(f1) > 1     
            return revLex([f1[i] for i in 1:length(f1)-1], [f2[i] for i in 1:length(f2)-1]);      
        else      
            return false;      
        end
   end
end;


E=[1,2,3,4,5,6,7,8,9];
Subs=subsets(Set(E),4);
Sub=[sort(x) for x in Subs];
PE=sort(Sub, lt=(x,y)->revLex(x,y));##Ordered subsets


NonHppCands=[];## List of the list of bases of matroids in nonHppcand file 
CandsData=[];##List of the matroid fingerprints
open("nonHPPCands.txt") do file
    for li in eachline(file)
        indexArray = findall( x -> x == "*", split(li,""))
	B=[PE[i] for i in indexArray]## List of bases of the matroid
        push!(NonHppCands,B)
        push!(CandsData,li)
    end
end;

using HomotopyContinuation

@polyvar x[1:9];##Polynomial ring

function BgP(B)##Function for finding the bases generating polynomial
    ##B is list of bases of the matroid
    r=length(B[1])
    h=0
    for i in 1:length(B)
        N=1
        b=B[i]
        for j in 1:r
            N=N*x[b[j]]
        end
        h=h+N
    end
    return h
end;

           
S1=[1,2,3,4,5,6,7,8];
S2=[1,2,3,4,5,6,7,8,9];
Jl=[(x,y) for x in S1, y in S2 if x!=y];
J=[]; ##List of indices (i,j)
for j in Jl
    if (j[2],j[1])∉ J
        push!(J,j)
    end
end;

##Creating files
pwd();
touch("NotHPP.txt");
touch("Undetected.txt");

for i in 1:length(NonHppCands)
    BM=NonHppCands[i]##List of bases of M
    DM=CandsData[i]##Fingerprint of M
    h=BgP(BM)##Bases generating polynomial of M
    for j in J
        dif1=differentiate([h],[x[j[1]]])
        dif2=differentiate([h],[x[j[2]]])
        dif12=differentiate(dif1,[x[j[2]]])
        rayl=dif1*dif2-h*dif12
        E=[1,2,3,4,5,6,7,8,9]
        Es=deleteat!(E,sort([j[1],j[2]]))
        
        for p in Es
            rayl=rayl+(1/1000)*x[p]^4
        end
        rayl=subs(rayl, x[Es[1]]=>1.0, x[j[1]]=>1.0, x[j[2]]=>1.0)
        f=[differentiate(rayl,[x[Es[2]]])[1],differentiate(rayl,[x[Es[3]]])[1],differentiate(rayl,[x[Es[4]]])[1],differentiate(rayl,[x[Es[5]]])[1],differentiate(rayl,[x[Es[6]]])[1],differentiate(rayl,[x[Es[7]]])[1]]
        
        result=solve(f)
        for k in real_solutions(result)
            kr=[1.0]
            for l in 1:length(k)
                lr=round(k[l]; digits=4)
                append!(kr,lr)
            end
            rayl=dif1*dif2-h*dif12
            ##rayl=subs(rayl, x[j[1]]=>1.0, x[j[2]]=>1.0)
            V=[x[i] for i in Es]
            F=subs(rayl[1],V=>kr)
            Ff=convert(Float64,F)
            if (Ff<0)##It found a negative point
                Of=open("NotHPP.txt","a")
                write(Of, "($DM,$kr,$j)\n");
                flush(Of)
                close(Of)
                @goto escape_label##Goes out of the loop 
            end

        end
    end
    oF=open("Undetected.txt","a")##When it couldn't find a negative point
    write(oF, "$DM\n");
    flush(oF)
    close(oF)
    @label escape_label
end;
    
    
    
    
    


