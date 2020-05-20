
include("./perov_functions.jl")
using Plots
using BenchmarkTools
using Distributions
using LaTeXStrings
using Measures

n = 2;
nk = 100;
system_original = make_symstem(n);
system = copy(system_original);

function get_g(de)
    n = 2
    system_original = make_symstem(n)
    system = copy(system_original)
    axis = [0, 1, 0.0]
    system[1, 2, 1] = rot_tetra_test_1(system[1, 2, 1], axis, pi)
    system[2, 1, 1] = rot_tetra_test_1(system[2, 1, 1], axis, pi)
    system[1, 2, 2] = rot_tetra_test_1(system[1, 2, 2], axis, pi)
    system[2, 1, 2] = rot_tetra_test_1(system[2, 1, 2], axis, pi)
    u_m = get_total_energy_mag(system, 0, 1)[1]
    d_m = get_total_energy_mag(system, 1, 0)[1]
    u_p = get_total_energy_mag(system_original, 0, 1)[1]
    d_p = get_total_energy_mag(system_original, 1, 0)[1]
    return (de - (d_p - d_m)) / (u_p - u_m)
end;

g = get_g(-0.095)
Ep =
    g * get_total_energy_mag(system, 0, 1)[1] +
    get_total_energy_mag(system, 1, 0)[1]
Em =
    g * get_total_energy_mag(system_original, 0, 1)[1] +
    get_total_energy_mag(system_original, 1, 0)[1];
println("g=$g  Check δE(0.095)=$(Ep-Em)") #= =0.095? =#

function run_sum(n, system, g1, g2)
    e = []
    a1 = []
    a2 = []
    a3 = []
    mag_val = []
    for i = 1:n
        axis = rand(Uniform(-1.0, 1.0), 3)
        theta = rand(Uniform(0, 2π))
        pos = rand(1:size(system)[1], 3)
        s, s1 = get_rotation_energy_1(system, g1, g2, pos, axis, theta)
        if s < 0
            system = s1
        end
        if i % 100 == 0
            vals = get_total_energy_mag(system, g1, g2)
            energy = vals[1]
            angles = vals[2] .^ 2
            push!(e, energy)
            push!(a1, angles[1])
            push!(a2, angles[2])
            push!(a3, angles[3])
            push!(mag_val, norm(vals[2]))
        end
    end
    return e, a1, a2, a3, mag_val, system
end;

g = get_g(-0.095)
n = 5;
nk = 100;
function animplots(i)
    gr()
    pl1=plot((1:i)*100,e[1:i],ylabel="energy",xlabel="iterations",color=RGB(226/255, 194/255, 144/255),label="",linewidth=2)
    pl2=plot((1:i)*100,a1[1:i],label=L"$s_x$",xlabel="iterations",ylabel="component",ylim=(-0.01,1.1),linewidth=1.5)
    plot!(pl2,(1:i)*100,a2[1:i],label=L"$s_y$",linewidth=1.5)
    plot!(pl2,(1:i)*100,a3[1:i],label=L"$s_z$",linewidth=1.5)
    pl3=plot((1:i)*100,mag_val[1:i],ylabel=L"$M$",xlabel="iterations",color=RGB(202/255, 255/255, 191/255),linewidth=2,label="",ylim=(-0.01,1.1))
    hline!(pl3,[1],linestyle=:dash,color="grey",label="")
    hline!(pl2,[1/3],linestyle=:dash,color="grey",label="")
    plot(pl1,pl3,pl2,size=(1200,400),dpi=100,bottom_margin=5mm,layout=@layout [a b c])
end

system = make_symstem(n);
e,a1,a2,a3,mag_val,system=run_sum(10000,system,1,0);
animplots(length(e))

function run_sum_at_T(n, system, g1, g2,T;capture=0)
    e = []
    a1 = []
    a2 = []
    a3 = []
    mag_val = []
    dir=[[1,0,0.],[0.,0,1],[1.,0,0]]
    lst=convert.(Int64,round.(LinRange(1,n,capture)))
    for i = 1:n
        axis = dir[rand(1:3)] #rand(Uniform(-1.0, 1.0), 3)
        theta = rand(Uniform(0, 2π))
        pos = rand(1:size(system)[1], 3)
        s, s1 = get_rotation_energy_1(system, g1, g2, pos, axis, theta)
        if s < 0 || rand() < exp(-s/T)
            system = s1
        end
        if capture >0
            if i in lst
                vals = get_total_energy_mag(system, g1, g2)
                energy = vals[1]
                angles = vals[2] .^ 2
                push!(e, energy)
                push!(a1, angles[1])
                push!(a2, angles[2])
                push!(a3, angles[3])
                push!(mag_val, norm(vals[2]))
            end;
        end;
    end
    vals = get_total_energy_mag(system, g1, g2)
    energy = vals[1]
    angles = vals[2] .^ 2
    push!(e, energy)
    push!(a1, angles[1])
    push!(a2, angles[2])
    push!(a3, angles[3])
    push!(mag_val, norm(vals[2]))
    if capture==0
        return e[1], a1[1], a2[1], a3[1], mag_val[1], system
    else
        return e, a1, a2, a3, mag_val, system
    end;
end;
system = make_symstem(n);
e = [];a1 = [];a2 = [];a3 = [];mag_val = []
T_range=LinRange(0,7e-2,20)
for T in T_range
    e_tmp,a1_tmp,a2_tmp,a3_tmp,mag_val_tmp,system_diff=run_sum_at_T(10000,system,1,0,T,capture=0);
    push!(e, e_tmp)
    push!(a1, a1_tmp)
    push!(a2, a2_tmp)
    push!(a3, a3_tmp)
    push!(mag_val, mag_val_tmp)
end;

pl_a=scatter(T_range,a1)
scatter!(pl_a,T_range,a2)
scatter!(pl_a,T_range,a3)
pl_e=scatter(T_range,e)
pl_m=scatter(T_range,mag_val)
plot(pl_e,pl_m,pl_a,size=(1200,400),dpi=100,bottom_margin=5mm,layout=@layout [a b c])
