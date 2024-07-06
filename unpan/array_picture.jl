begin
    using LinearAlgebra
    using Plots
    using JLD2
    using PyCall
    using Printf
end


begin
    #入力してね
    diry = "fixed2α"
    Npick = 1
end

begin
    @load "./$diry/Eovl.$diry.jld2" Eout ovlout proHout bareH bareovl
    #@load "./$diry/coor$diry.jld2" coord gate    
end

R = zeros(Float64,12)
energy = zeros(Float64,12)
for i in 1:12
    R[i] = i 
    energy[i] = Eout[i][1]
end

output ="./2a.png"
p = plot([R],[energy], label="", title="α-α distance v.s. <E>",xlabel="α-α distance(fm)", ylabel="energy(MeV)")
savefig(p,output)