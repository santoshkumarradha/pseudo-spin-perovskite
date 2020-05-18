#---- Imports
using LinearAlgebra
using StaticArrays
using Pkg
using PyCall
using Plots
Plots.plotlyjs()
using Distances
using Statistics
theme(:dark)
p = pyimport("pymatgen")
using Distributions


#= Defining structure for tetra
=#
struct tetra
    ge::Array{Float64,1}
    cl1::Array{Float64,1}
    cl2::Array{Float64,1}
    cl3::Array{Float64,1}
end


#check equavality of two struc 
function check_equal(t1::tetra,t2::tetra)
    if t1.ge==t2.ge && t1.cl1==t2.cl1  && t1.cl2==t2.cl2 && t1.cl3==t2.cl3
        return true
    else
        return false
        end;
    end;

#= Mat mul for rotation
=#
icross(b) = copy(Transpose(hcat([cross(Matrix(1.0I, 3, 3)[:,i],b) for i in 1:3]...)));
anchor(tetra) = (tetra.ge+tetra.cl1+tetra.cl2+tetra.cl3)/4

function Mdot(a1,a2)
    #=redefing dot product like numpy for matrix=#
    a1_1=copy(Transpose(a1))'
    return [dot(a1_1[i,:],a2) for i in 1:3]
end


function rot(coords,anc,axis,theta)
    #=rotate coords by theta wrt anc with axis as axis=#
    theta %= 2 * pi
    rm=exp(icross(axis/norm(axis))*theta)
    val=Mdot(rm,(coords-anc))+anc
    return val
    end;


#ROtating a tetra
function rot_tetra_test_1(t::tetra,axis::Array{Float64,1},theta)
    #=rotate all atoms in tetra by theta with axis = axis=#
    anc=anchor(t)
    ge=rot(t.ge,anc,axis,theta)
    cl1=rot(t.cl1,anc,axis,theta)
    cl2=rot(t.cl2,anc,axis,theta)
    cl3=rot(t.cl3,anc,axis,theta)
    return tetra(ge,cl1,cl2,cl3)
    end;

#---PBC index
function pbc_1(system,index)
    eval_1(i,n)= ((abs(i)-1) ÷ n)*((sign(i)+1) ÷ 2) - ((sign(i-1)-1) ÷ 2)*((i-n) ÷ n)
    index_1(i,n)=((abs(i+n-1) % n)+1)
    unit_size=size(system)[1]
    a=(-system[1,1,1].ge+system[1,1,2].ge)[3]*unit_size
    uc=eval_1.(index,unit_size)*a
    index_mod=index_1.(index,unit_size)#index .% (unit_size+1)
    return index_mod,uc
    end;

#----get the tetra of a system at given index with PBC
function system_at_index(system,index)
    index_mod,uc=pbc_1(system,index)
    tetra_tmp=system[index_mod[1],index_mod[2],index_mod[3]]
    ge=tetra_tmp.ge+uc
    cl1=tetra_tmp.cl1+uc
    cl2=tetra_tmp.cl2+uc
    cl3=tetra_tmp.cl3+uc
    return tetra(ge,cl1,cl2,cl3)
    end;




# NN algorithm
function get_nn_2(system,pos)
    tmp=Array{Float64, 1}[]
    for i in -1:1
        for j in -1:1
            for k in -1:1
                index_tmp=pos+[i,j,k]
                temp=system_at_index(system,index_tmp)
                push!(tmp,temp.cl1)
                push!(tmp,temp.cl2)
                push!(tmp,temp.cl3)
            end
        end
    end
    return tmp
end

# Mean and Var from NN
function get_mean_var(system,pos,return_type="var")
    # Get the mean or variance of a position of lattice  wrt bond distance#
    sys=system[pos[1],pos[2],pos[3]]
    distance_1=colwise(Euclidean(), sys.ge, copy(hcat(get_nn_2(system,pos)...)))
    nn=sort(distance_1)[1:6]
    var_sys=var(nn)
    mean_sys=mean(nn)
    if return_type=="mean"
        return mean_sys
    else
        return var_sys
    end
end


# Make the system from CsSiI2.cif and pymatgen

function make_symstem(n=2)
    cssii2=p.Structure.from_file("Cssii2.cif")
    struc1=cssii2.copy()
    struc1.make_supercell([[n,0,0],[0,n,0],[0,0,n]])
    dist=cssii2.get_distance(1,3)+.1
    positions=Array{Int64, 1}[]
    for i in struc1
        if i.species_string=="Si"
            push!(positions,struc1.get_neighbor_list(dist,[i])[2])
        end
    end
    system=Array{tetra,3}(undef,n,n,n);
    cnt=1
    for i in 1:n
        for j in 1:n
            for k in 1:n
                ge=get(struc1,reverse(positions[cnt])[1]).coords
                cl1=get(struc1,reverse(positions[cnt])[2]).coords
                cl2=get(struc1,reverse(positions[cnt])[3]).coords
                cl3=get(struc1,reverse(positions[cnt])[4]).coords
                system[i,j,k]=tetra(ge,cl1,cl2,cl3)    
                cnt+=1
            end
        end
    end
    return system
end

# Function to get dipole orientation of a tetrahedron

function get_dipol_vec(tetra::tetra)
    return (tetra.ge - (tetra.cl1+tetra.cl2+tetra.cl3)/3)/norm(tetra.ge - (tetra.cl1+tetra.cl2+tetra.cl3)/3)
end

# Get the distance  vector between two tetrahedrons supplying pl will plot it in pl
function get_dipole_dist(tetra1::tetra,tetra2::tetra,pl=nothing)
    r1=(tetra1.ge+tetra1.cl1+tetra1.cl2+tetra1.cl3)/4
    r2=(tetra2.ge+tetra2.cl1+tetra2.cl2+tetra2.cl3)/4
    if pl == nothing
        return r1-r2
    else
        plot_tetra(tetra1,pl)
        plot_tetra(tetra2,pl)
        plot!(pl,[i[1] for i in [r1,r2]],[i[2] for i in [r1,r2]],[i[3] for i in [r1,r2]],color="red",linewidth=3,label="connecting vector r12")
        end;
    end;



#=--- dipole energy between two tetra hedrons
BigFloat is needed for accuracy. Been struggling for ages without that ! sigh. =#

# function get_dipole_energy_between_tetra(tetra1::tetra,tetra2::tetra)
#     r12=[BigFloat(i) for i in get_dipole_dist(tetra1,tetra2)]*1.0E-10
#     d1=[BigFloat(i) for i in get_dipol_vec(tetra1)]*1.0E-10
#     d2=[BigFloat(i) for i in get_dipol_vec(tetra2)]*1.0E-10
#     four_pe0= 9 * 10^9 #1/4πe0 in N⋅m^2⋅C^−2
#     return BigFloat(( dot(d1,d2) - 3*(dot(d1,r12)*dot(d2,r12))/norm(r12)^2 ) * four_pe0 * norm(r12)^-3)
#     end;


#This seems to be better ? 
function get_dipole_energy_between_tetra(tetra1::tetra,tetra2::tetra)
    r12=get_dipole_dist(tetra1,tetra2)
    r12=r12
    d1=get_dipol_vec(tetra1)
    d2=get_dipol_vec(tetra2)
    four_pe0= 9 * 10^9 #1/4πe0 in N⋅m^2⋅C^−2
    #return BigFloat(( dot(d1,d2)*norm(r12)^-2 - 3*(dot(d1,r12)*dot(d2,r12))*norm(r12)^-4 ))
    return -dot(d1,d2)+3*(dot(d1,r12)*dot(d2,r12))
    end;

#= Get dipole energy of a given system at pos
Need to check the other algorithm giving weired results=#
function get_dipole_energy(system,pos)
    nn=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    d=0
    tetra1=system[pos[1],pos[2],pos[3]]
    for i in nn
        pos_new=pos+i
        temp=system_at_index(system,pos_new)
        d+=get_dipole_energy_between_tetra(tetra1,temp)
        end;
    return d
    end;

#Rotate the system
function rotate_system(system,pos,axis,theta)
     system[pos[1],pos[2],pos[3]]=
            rot_tetra_test_1(system_at_index(system,pos),axis,theta);
    return system
    end;


#= Rotation total energy cost calculation
=#
function get_rotation_energy_1!(system,g1,g2,pos,axis,theta)
    system_test=copy(system)
    system_test[pos[1],pos[2],pos[3]]=
            rot_tetra_test_1(system_at_index(system,pos),axis,theta);
    nn=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    u=0
    for i in nn
        ind,_=pbc_1(system,pos+i)
        print(ind)
        u+=(get_mean_var(system_test,ind,"var")-get_mean_var(system,ind,"var"))
    end
    d=get_dipole_energy(system_test,pos)-get_dipole_energy(system,pos)
    system[pos[1],pos[2],pos[3]]=system_test[pos[1],pos[2],pos[3]]
    return g1*u+g2*d
    end;
# this function returns the system. Which is helpful for updating the value
function get_rotation_energy_1(system,g1,g2,pos,axis,theta)
    system_test=copy(system)
    system_test[pos[1],pos[2],pos[3]]=
            rot_tetra_test_1(system_at_index(system,pos),axis,theta);
    nn=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    u=0
    for i in nn
        ind,_=pbc_1(system,pos+i)
        u+=(get_mean_var(system_test,ind,"var")-get_mean_var(system,ind,"var"))
    end
    d=get_dipole_energy(system_test,pos)-get_dipole_energy(system,pos)
    return g1*u+g2*d,system_test
    end;


#= Total energy and total Magnetic direction (magnatization)
of the system normalized over total tetrahedrons (it seems costly)=#
function get_total_energy_mag(system,g1,g2)
    system_size=size(system)[1]
    d=0
    u=0
    m=[0,0,0]
    for i in 1:system_size
        for j in 1:system_size
            for k in 1:system_size
                pos=[i,j,k]
                d+=get_dipole_energy(system,pos)
                u+=get_mean_var(system,pos,"var")
                m+=get_dipol_vec(system[i,j,k])
                end;
            end;
        end;
    u/=prod(size(system))
    d/=prod(size(system))
    m/=prod(size(system))
    return g1*u+g2*d,m
end











#=------------------
Plotting Tools
-------------------=#
function plot_tetra(tetra::tetra,pl)
    pos1=[tetra.ge,tetra.cl1,tetra.cl2,tetra.cl3]
    scatter!(pl,[pos1[i][1] for i in [2,3,4]],
        [pos1[i][2] for i in [2,3,4]],
        [pos1[i][3] for i in [2,3,4]],color="darkblue",label="")
    scatter!(pl,[pos1[i][1] for i in [1]],
        [pos1[i][2] for i in [1]],
        [pos1[i][3] for i in [1]],color="darkred",label="")
    for i in 2:4
        plot!(pl,[pos1[1][1],pos1[i][1]],
            [pos1[1][2],pos1[i][2]],
            [pos1[1][3],pos1[i][3]],linewidth=5,color="white",label="")
    plot!(pl,[pos1[i][1] for i in [2,3,4,2]],
            [pos1[i][2] for i in [2,3,4,2]],
            [pos1[i][3] for i in [2,3,4,2]],linewidth=5,
            color="grey",linestyle=:dash,label="")
    end
end









#=---- example
n=4
system=make_symstem(n);
pl1=plot(ticks = false)
system_test=copy(system)
pos_1=[rand(1:4),rand(1:4),1]
    axis=[0,0,1.]
    theta=pi
    system_test[pos_1[1],pos_1[2],pos_1[3]]=
            rot_tetra_test_1(system[pos_1[1],pos_1[2],pos_1[3]],axis,theta);

for i in system_test
    plot_tetra(i,pl1)
end
plot!(pl1,size=(900,900))
display(pl1)
=#


