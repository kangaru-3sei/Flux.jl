#8heHSPfull配位における図の作成

begin
    using LinearAlgebra
    using Plots
    using JLD2
    #using PyCall
    using Printf
end


#include("C:/Users/chikk/OneDrive - Kyoto University/analysis/GCM_analysis.jl")
cd("unpan")
begin
    #入力してね
    diry = "8he.ver6"
    Npick = 5
    @load "./$diry/Eovl.$diry.jld2" proH Etrn wf H ol trn evol norvec tN Hvec_trn 
    @load "./$diry/coor$diry.jld2" HR gate
    N = size(H,1)
    coord = HR
    include("GCM_analysis_func.jl")
    c = GCM_keisu(evol,Hvec_trn,norvec)
end

#closure配位とのオーバーラップ
abs2(wf[1,1])
abs2(wf[1,2])
8.6/48.3
(29.1-16.9)/29.1
 0.4
 0.7444886258389696
 0.6325817464895113
 0.8647976791292165
Etrn[1]-Etrn[2]
#ls0を抜き出す
begin
    begin
        #入力してね
        diry = "12CM6.ls0"
        Npick = 5
    end
    begin
        @load "./$diry/Eovl.$diry.jld2" proH Etrn wf H ol trn evol norvec tN Hvec_trn 
        @load "./$diry/coor$diry.jld2" HR gate    
    end
    proHls0 = proH
    Etrnls0 = Etrn 
    wfls0= wf
    Hls0=H
    olls0=ol
    evolls0=evol
    norvecls0=norvec
    Hvec_trnls0=Hvec_trn
end
=#

Rname = "12CHSP.rad"
begin
    include("./GCM_analysis.radii.jl")
    @load "./$Rname/radii.$Rname.jld2" radius
    rad = radius   
end

@load "./$diry/radii.$diry.rad.jld2" radius
rad = radius
Nrad - Prad
2*rad - Prad -  Nrad
Prad - Nrad
#@save "./$diry/radii.1bodyop.jld2" rad Nrad Prad
8 .* rad - 2 .* Prad - 6 .* Nrad
#ある物理量の期待値を求める。状態は固有所帯iの fullのみ
#radiusは求めたい物理量をgcm配位で挟んだもの
expecval_full(rad,2)

7.680145*12
#３体の配置のみを抜き出す
begin
    coordp = HR[gate,:]
    wfp = real.(wf[gate,:])
end

begin
    sc_interpolate = pyimport("scipy.interpolate")
    np = pyimport("numpy")
end

#θ=90のみを記録

#=
A = [1 2 3; 4 5 6; 7 8 9]  # 3x3行列
mask = A[:, 1] .> 5  # 第1列の値が5より大きいかどうかを示すブール配列
B = A[mask, :]  # 条件を満たす行だけを抽出
=#

#上を参照に抜き出す

maskθ = coordp[:,3] .== π/2
maskY0 = coordp[:,2] .== 0
mask = maskθ .| maskY0

#=
scatter([Xeqd],[Yeqd],yrange=(0,12),aspect_ratio=:equal)
=#
#Y,Xの座標でプロット
Y = coordp[:,1] .* sin.(coordp[:,2])
X = coordp[:,1] .* cos.(coordp[:,2]) 

Xeqd = X[mask][:]
Yeqd = Y[mask][:]
wfeqd = wfp[mask,:]
#いつもの
for i in 1:Npick    # ith state
    wfp = abs2.(wfeqd[:,i])
    
    # 補間を行う新しいグリッドの生成
    lingrid_y = LinRange(minimum(Yeqd), maximum(Yeqd), 100)
    lingrid_x = LinRange(minimum(Xeqd), maximum(Xeqd), 100)
    grid_x, grid_y = np.meshgrid(lingrid_x, lingrid_y)
    
    # griddataを使用して補間し、範囲外の点は0で埋める
    zi = sc_interpolate.griddata((Xeqd, Yeqd), wfp, (grid_x, grid_y), method="cubic", fill_value=0)

    output = joinpath(diry, "wfXYfull.st$i.png")
    p = contourf(lingrid_x, lingrid_y, zi, xlabel="X", ylabel="Y", xrange=(0,12), yrange=(0,12), aspect_ratio=:equal, c=:thermal, title="", tickfontsize=15, labelfontsize=15)
    savefig(p, output)
end

    #実空間のほう
    #=
    lingrid_r2 = LinRange(minimum(r2eqd), maximum(r2eqd), 100)
    grid_x, grid_r2 = np.meshgrid(lingrid_x, lingrid_r2)
    zi = sc_interpolate.griddata((Xeqd, r2eqd), wfp, (grid_x, grid_r2), method="cubic", fill_value=0)
    output = "wfR.st$i.r2.png"
    p = contourf(lingrid_x, lingrid_r2, zi, xlabel="R1", ylabel="R2",  aspect_ratio=:equal,xrange=(0,12/√2), yrange=(0,12), c = :thermal, title="wave function of $i th 0⁺ state")
    savefig(p, output)
    =#


#6he+2n配位の場合
maskshell = gate[:] .== false
maskshell[1] = false
co6he = coord[maskshell,:]
wf6he = wf[maskshell,:]
R6he = co6he[:,1]
#同じRの中で１番目、２番目の配位を選び出す
fsgate1 = fill(true, length(R6he))
for i in 1:length(R6he)
    if i%2 == 0
        fsgate1[i] = false
    end
end
R6hep = R6he[fsgate1,:]     #plotに使用
fsgate2 = fsgate1[:] .== false
wf6he1 = wf6he[fsgate1,:]   #一つ目の回転状態
wf6he2 = wf6he[fsgate2,:]   #二つ目の〃

#各状態に対する2配位の波動関数
for i in 1:Npick
    if i == 1 || i == 2
    wfp1 = real.(wf6he1[:,i])
    wfp2 = real.(wf6he2[:,i])

    output = "./$diry/wf6he.st$i.png"
    p = plot([R6hep], [wfp1], type=:scatter, markershapes=:circle, framestyle=:origin, color=:blue, ylims=(-0.8,0.8))
    plot!([R6hep], [wfp2], markershapes=:circle, color=:red)
    savefig(p, output)
    end
end

#配位1における状態とのオーバーラップ
begin
    output = "./$diry/wf6he1.png"
    output2 = "./$diry/wf6he2.png"
    p = plot(framestyle=:origin, ylims=(-0.9,0.9))
    q = plot(framestyle=:origin, ylims=(-0.9,0.9))
    for st in 1:Npick
        if st == 1 || st == 2
            wfp = real.(wf6he1[:,st])
            plot!(p, R6hep, wfp, markershapes=:circle, line=:solid)
            wfp = real.(wf6he2[:,st])
            plot!(q, R6hep, wfp, markershapes=:circle, line=:solid)
        end
    end
    savefig(p, output)
    savefig(q,output2)
end


#closure, 6he+2n, θ90 の配位のゲートを作成
begin
    #closure配位のgate
    begin
        clgate = fill(false, size(H,1))
        clgate[1] = true
        clgate = BitVector(clgate)
    end

    #6he+2n配位
    begin
        gate6he = .!gate
        gate6he[1] = false
    end  

    #lsが入ったもののみ
    begin
        gatels = clgate .| gate6he
    end

    #θ=90の配位
    begin
        mask90 = coord[:,3] .== π/2
        
        cogate = coord[:,2] .== 0    #α=0かつgate==trueのものが3体が直線に並ぶ配位
        lingate = gate .& cogate

        gate90 = mask90 .| lingate      #θ=90配位のgate
    end
end

gate
#subspaceとのオーバーラップ。subspaceでのnatural basisとj th state in full spaceとのオーバーラップの値を出力 
#inputはgate:subspace指定 のみ。
calcu_subovl(gate)
calcu_subovl(gate6he)
gatemicro,count = cluster_distance_cutoff_8he(0.2)
calcu_subovl(gatemicro)

#gateから指定したRのみをさらに絞り込むゲートを返す
function gateRtrncate(gate::BitVector, coord, R)  
    gateR = fill(false, length(gate))
    for i in axes(coord,1)
        if coord[i,1] == R
            gateR[i] = true
        end
    end
    gateRout = gateR .& gate
    return gateRout
end

#時間がかかるので注意！
#6he subspace とのovlp
resall = Vector{Vector{Float64}}()
for R in 1:12
    gateR6he = gateRtrncate(gate6he, coord, R)
    res = calcu_subovl(gateR6he)
    push!(resall, res)
end

#3body 配位とのsqared ovlps
resall = Vector{Vector{Float64}}()
for R in 1:12
    gate3R = gateRtrncate(gate, coord, R)
    res = calcu_subovl(gate3R)
    push!(resall, res)
end

#そのプロット。collective spaceで対角化した各状態とのovlpと、その和をとったものをそれぞれ出力。
for i in 1:2    #i th 固有状態　of Ψfull
    p = plot(yrange=(0,1.0), xticks=(0.0:1.0:13.0),title="sqrd. ovlp. of st$i with 3-body subspace", xlabel="R (fm)")
    output = joinpath(diry, "ovlp.3body.st$(i).png")
    line1 = Vector{Float64}(undef, 12)
    yoko = Vector{Int8}(undef, 12)
    
    for j in 1:12
        res = resall[j][i]
        line1[j] = res
        yoko[j] = j
    end
    plot!(p, yoko, line1, type=:scatter, markershapes=:circle, label="",legendfontsize = 15, tickfontsize = 15)
    # 各点の隣にその値を書く
    for j in 1:12
        offset = 0.05 # オフセットの値を調整する
        annotate!(p, yoko[j], line1[j] + offset, text(string(round(line1[j], digits=3)), 8))
    end
    savefig(p,output)
    display(p)

    #Hで対角化した各状態とのovlpの和をとる
    if i==1 #出力の絵を１枚にする。
        global q = plot(yrange=(0,1.0),xticks=(0.0:1.0:13.0),title="sqrd. ovlp. with 3-body subspace", xlabel="R (fm)")
        global r = plot(yrange=(0,1.0),xticks=(0.0:1.0:13.0),title="sqrd. ovlp. with 3-body subspace", xlabel="R (fm)")
        labelname = "g.s."
    else 
        labelname = "2nd"
    end
    outputgentle = joinpath(diry, "ovlp.3body.yasashi.png")
    output = joinpath(diry, "ovlp.3body.png")
    plot!(q, yoko, line1, type=:scatter, markershapes=:circle, label=labelname, legendfontsize = 15, tickfontsize = 15)
    plot!(r, yoko, line1, type=:scatter, markershapes=:circle, label=labelname, legendfontsize = 15, tickfontsize = 15, labelfontsize=15)
    savefig(r, output)
    for j in 1:12
        offset = 0.05 # オフセットの値を調整する
        annotate!(q, yoko[j], line1[j] + offset, text(string(round(line1[j], digits=3)), 8))
    end
    savefig(q,outputgentle)
    display(q)
end

function matzero(mat::Hermitian{ComplexF64, Matrix{ComplexF64}},mask)
    #matとgateのサイズが違う場合、警告を鳴らす
    if length(mask) != size(mat,1)
        println("gateと行列のサイズが違いますわ...")
        return 0
    end

    # gateのfalseの位置を見つける
    line_remove = findall(.!mask)

    mat = Matrix(mat)   #hermiteのままでは0にできない。

    # 指定した列を0に設定
    mat[:, line_remove] .= 0
    # 指定した行を0に設定
    mat[line_remove, :] .= 0

    #hermiteにする
    mat = Hermitian(mat)
    return mat
end


#subspaceで対角化してえられた固有状態と特定の配位とのovlp
#inputはsubspaceを形成するgateoと、ovlpをとる配位を指定するi
function ovlp_wconf_osub(bareH::Hermitian{ComplexF64, Matrix{ComplexF64}},bareol::Hermitian{ComplexF64, Matrix{ComplexF64}}, gateo, i::Int, trn)
    n = size(bareH,1)   #全空間の次元
    #gateでsubspaceの構成
    Ho = matzero(bareH, gateo)
    olo = matzero(bareol, gateo)
    Hovec, eigval_olo, olovec = Hvec_trn_sub(Ho, olo, trn)
    subdim = size(Hovec,1)  #全固有状態について一応計算

    Hovec = Hvec_kakucho(Hovec, n)  #n*nにHovecを拡張
    start = n-subdim+1  #Hovecが始まるところ
    resvec = zeros(Float64,subdim)
    for st in start:n 
        res = 0.0 +0.0im
        for k1 in start:n
            atama = Hovec[k1,st]/sqrt(eigval_olo[k1])
            siri = 0.0+0.0im
            for a in 1:n
                if gateo[a] == true
                    siri += olovec[a,k1]*bareol[i,a]
                end
            end
            res += atama*siri
        end
        resvec[st-start+1] = abs2(res)/bareol[i,i]
    end
    totres = sum(resvec)
    return totres, resvec
end

ovlp_osub(gate,i) = ovlp_wconf_osub(H,ol,gate,i,trn)
gatefull = trues(490)
#closure配位とのovlp
ovlpclfull = ovlp_osub(gate,1)
gatewols = .!(gate6he .| clgate)
ovlpclNols = ovlp_osub(gatewols,1)
ovlpcl6he = ovlp_osub(gate6he,1)
gatels = gate6he .| clgate
ovlpwls = ovlp_wconf_osub(H,ol,gatels,1,1.0e-18)
gate31 = .!gate6he
ovlpw31 = ovlp_wconf_osub(H,ol,gate31,1,trn)
gatels90 = gate90 .| clgate .| gate6he
ovlpls90 = ovlp_osub(gatels90,1)
gatels90 = gate90 .| clgate
ovlpls90 = ovlp_wconf_osub(H,ol,gatels90,1,trn)
gate23 = gate6he .| gatewols
ovlp23 = ovlp_osub(gate23,1)



#↑を拡張してsubspaceとの overlapを求める
#例えばls90と3-bodyとのovlpとか。
function ovlp_wall_osub(bareH::Hermitian{ComplexF64, Matrix{ComplexF64}},bareol::Hermitian{ComplexF64, Matrix{ComplexF64}}, gateo::BitVector, trn)
    n = size(bareH,1)   #全空間の次元
    #gateでsubspaceの構成
    Ho = matzero(bareH, gateo)
    olo = matzero(bareol, gateo)
    Hovec, eigval_olo, olovec = Hvec_trn_sub(Ho, olo, trn)
    subdim = size(Hovec,1)

    Hovec = Hvec_kakucho(Hovec, n)  #n*nにHovecを拡張
    start = n-subdim+1  #Hovecが始まるところ
    stop = start + 3    #手間を省くため、4th状態までに限定

    #作成したsubspaceの固有状態と全配位とのovlp
    resvec = zeros(Complex{Float64},n,4)  
    
    for b in 1:n 
        for st in start:stop
            res = 0.0 +0.0im
            for k1 in start:n
                atama = Hovec[k1,st]/sqrt(eigval_olo[k1])
                siri = 0.0+0.0im
                for a in 1:n
                    if gateo[a] == true
                        siri += olovec[a,k1]*bareol[b,a]
                    end
                end
                res += atama*siri
            end
            resvec[b,st-start+1] = res
        end
    end
    return resvec
end


#with subspaceのほうの、gcm係数を求める
#↑をのループの位置を入れ替えただけ
function wsub_keisu(bareH::Hermitian{ComplexF64, Matrix{ComplexF64}},bareol::Hermitian{ComplexF64, Matrix{ComplexF64}}, gatew::BitVector, trn)
    #gatewでwsubspaceの構成
    Hw = matremove(bareH, gatew)
    olw = matremove(bareol, gatew)
    n = size(Hw,1)   #sub空間の次元
    Hwvec, eigval_olw, olwvec = Hvec_trn_sub(Hw, olw, trn)
    subdim = size(Hwvec,1)

    elimstart = length(eigval_olw) - subdim + 1
    eigval_olw = eigval_olw[elimstart:end]

    #wsubspaceのgcm係数
    resvec = zeros(Complex{Float64},n,subdim)  

    for st in 1:subdim  #こちらはすべての固有状態について足す。
        for a in 1:n
            res = 0.0 +0.0im
            for k1 in 1:subdim
                res += Hwvec[k1,st]/sqrt(eigval_olw[k1])*olwvec[a,k1]
            end
            resvec[a,st] = res
        end
    end    
    return resvec
end

function ovlp_wsub_osub(bareH::Hermitian{ComplexF64, Matrix{ComplexF64}},bareol::Hermitian{ComplexF64, Matrix{ComplexF64}}, gatew::BitVector, gateo::BitVector, trn)
    keisu = wsub_keisu(bareH, bareol, gatew, trn)
    ovlp_wconf =  ovlp_wall_osub(bareH, bareol, gateo, trn)
    
    n = size(bareH,1)
    dimw = size(keisu,2)
    resvec = zeros(Float64,4)
    for st in 1:4   #面倒なので4状態のみ
        res = 0.0+0.0im
        for i in 1:dimw
            a_taio = 1
            for a in 1:n
                if gatew[a] == true
                    res += keisu[a_taio,i]*ovlp_wconf[a,st]
                    a_taio += 1
                end
            end
        end
        resvec[st] = abs2(res)
    end
    return resvec
end

#r.m.s. radiiなどの期待値
evol
Hvec_trn
norvec
radius
c = zeros(Complex{Float64},490,tN)
for j in 1:490
    for i in 1:tN
        res = sum(Hvec_trn[k,i]/√evol[k]*norvec[j,k] for k in 1:tN)
        c[j,i] = res
    end
end
#test: overlapはOK
try
    ovlp = zeros(Float64, 10,10)
for i in 1:10
    for j in 1:10
        ovlp[i,j] = abs2(sum(conj(c[k,j])*c[l,i]*ol[k,l] for k in 1:N, l in 1:N))
        if i==j
            if abs(1-ovlp[i,j]) > 0.001
                println(i)
            end
        else
            if ovlp[i,j] > 0.001
                println(i, ", ",j)
            end
        end
    end
end 
catch
end

#r.m.s. radii of ith state
res = 0.0+0.0im
i=1
N=size(H,1)
res = sum(conj(c[k,i])*c[l,i]*rad[k,l] for k in 1:N, l in 1:N)
res = sqrt(res)

0.235*2*197^2/940
#遷移強度
i = 1
j = 2
c
res = sum(conj(c[k,j])*c[l,i]*radius[k,l] for k in 1:N, l in 1:N)
res = abs(res) * 12

#ある物理量の期待値を求める。状態は固有所帯iの fullのみ
#radiusは求めたい物理量をgcm配位で挟んだもの
expecval_full(rad,1)^2
Etrn

keisu = wsub_keisu(H, ol, gate6he, trn)
ovlp_wconf =  ovlp_wall_osub(H, ol, gate, trn)
ovlp6he = ovlp_wsub_osub(H,ol, gate6he, gate6he, trn)

#２クラスター間距離が特定の値より小さいsubspaceを作成するgate
    #8heの場合
    coord[:,1]
function cluster_distance_cutoff_8he(r,Ntot=N,gate3=gate)
    gatekyori = trues(Ntot)
    count = length(gate3) + 24
    for i in eachindex(gate3)
        if gate3[i] == true
            R = coord[i,1]
            α = coord[i,2]
            θ = coord[i,3]
            R1 = R*cos(α)
            R2 = R*sin(α)/√2
            x2 = -0.5*R2*sin(θ)
            z2 = 0.5*(R1-R2*cos(θ))
            xα = -x2
            zα = 0.50*R2*cos(θ)
            kyori_α2n = sqrt((zα-z2)^2+(xα-x2)^2)
            if kyori_α2n >= r && R1 >= r
                gatekyori[i] = false
                count -= 1
            end
        else 
            gatekyori[i] = false
            count -= 1
        end

        #6he+2n配位
            gate6he = .!gate
            gate6he[1] = false
            #gatekyori = gatekyori .| gate6he
    end
    return gatekyori, count
end
    #12Cの場合
    function cluster_distance_cutoff_12C(r,Ntot=N,gate3=gate,zahyo=coord)
        gatekyori = trues(Ntot)
        count = length(gate3)
        for i in eachindex(gate3)
            if gate3[i] == true
                R = coord[i,1]
                α = coord[i,2]
                θ = coord[i,3]
                R1 = R*cos(α)/√2
                R2 = R*sin(α)*sqrt(3/8)
                x1 = 2/3*R2*sin(θ)
                x = -1/3*R2*sin(θ)
                z1 = 2/3*R2*cos(θ)
                z2 = (0.50*R1-1/3*R2*cos(θ))
                kyori_2α = sqrt((x-x1)^2+(z1-z2)^2)
                if kyori_2α ≥ r && R1 ≥ r
                    gatekyori[i] = false
                    count = count - 1
                end
            else 
                gatekyori[i] = false
                count = count - 1
            end
        end
        return gatekyori, count
    end    
sub2dis, count = cluster_distance_cutoff_8he(6)
calcu_subovl(sub2dis)
resall = Vector{Vector{Float64}}()
for R in 1:12
    gateRsub = gateRtrncate(sub2dis, coord, R)
    res = calcu_subovl(gateRsub)
    push!(resall, res)
end

for i in 1:2    #i th 固有状態　of Ψfull
    p = plot(yrange=(0,1.0), xticks=(0.0:1.0:13.0),title="sqrd. ovlp. of st$i with subspace cut by 2α distance", xlabel="R (fm)")
    output = joinpath(diry, "ovlp.2cluster.distance.st$(i).png")
    line1 = Vector{Float64}(undef, 12)
    yoko = Vector{Int8}(undef, 12)
    
    for j in 1:12
        res = resall[j][i]
        line1[j] = res
        yoko[j] = j
    end
    plot!(p, yoko, line1, type=:scatter, markershapes=:circle, label="",legendfontsize = 15, tickfontsize = 15)
    # 各点の隣にその値を書く
    for j in 1:12
        offset = 0.05 # オフセットの値を調整する
        annotate!(p, yoko[j], line1[j] + offset, text(string(round(line1[j], digits=3)), 8))
    end
    savefig(p,output)
    display(p)

    #Hで対角化した各状態とのovlpの和をとる
    if i==1 #出力の絵を１枚にする。
        global q = plot(yrange=(0,1.0),xticks=(0.0:1.0:13.0),title="sqrd. ovlp. with subspace cut by α-2n distance", xlabel="R (fm)")
        global r = plot(yrange=(0,1.0),xticks=(0.0:1.0:13.0),title="sqrd. ovlp. with subspace cut by α-2n distance", xlabel="R (fm)")
        labelname = "g.s."
    else 
        labelname = "2nd"
    end
    outputgentle = joinpath(diry, "ovlp.2cluster.distance.yasashi.png")
    output = joinpath(diry, "ovlp.2cluster.distance.png")
    plot!(q, yoko, line1, type=:scatter, markershapes=:circle, label=labelname, legendfontsize = 15, tickfontsize = 15)
    plot!(r, yoko, line1, type=:scatter, markershapes=:circle, label=labelname, legendfontsize = 15, tickfontsize = 15, labelfontsize=15)
    savefig(r, output)
    for j in 1:12
        offset = 0.05 # オフセットの値を調整する
        annotate!(q, yoko[j], line1[j] + offset, text(string(round(line1[j], digits=3)), 8))
    end
    savefig(q,outputgentle)
    display(q)
end

#クラスター距離に応じたサブスペースとのオーバーラップを見る
#12Cの場合
resall = Vector{Vector{Float64}}()
darray = Float64[]
dmax = 6
dmin = 0.05
d = dmin
while d ≤ dmax 
    dgate, count = cluster_distance_cutoff_8he(d)
    resvec = calcu_subovl(dgate)
    push!(resall, resvec)
    push!(darray, d)
    if d ≤ 2
        d = d + 0.1
    else
        d = d + 0.5
    end
end
resall[20]
for i in 1:2    #i th 固有状態　of Ψfull
    p = plot(yrange=(0,1.0), xticks=(0.0:1.0:13.0),title="sqrd. ovlp. of st$i with subspace cut by 2n-2n distance", xlabel="dineutron distance (fm)")
    output = joinpath(diry, "ovlp.2cluster.distance.subspacest$(i).png")
    line1 = Float64[]
    yoko = darray
    
    for j in axes(resall,1)
        res = resall[j][i]
        push!(line1, res)
    end
    plot!(p, yoko, line1, type=:scatter, markershapes=:circle, label="",legendfontsize = 15, tickfontsize = 15)
    # 各点の隣にその値を書く
    for j in axes(resall,1)
        offset = 0.05 # オフセットの値を調整する
        annotate!(p, yoko[j], line1[j] + offset, text(string(round(line1[j], digits=3)), 8))
    end
    savefig(p,output)
    display(p)

    #Hで対角化した各状態とのovlpの和をとる
    if i==1 #出力の絵を１枚にする。
        global q = plot(yrange=(0,1.0),xticks=(0.0:1.0:13.0),title="sqrd. ovlp. with subspace cut by 2n-2n distance", xlabel="dineutron distance (fm)")
        global r = plot(yrange=(0,1.0),xticks=(0.0:1.0:13.0),title="sqrd. ovlp. with subspace cut by 2n-2n distance", xlabel="dineutron distance (fm)")
        labelname = "g.s."
    else 
        labelname = "2nd"
    end
    outputgentle = joinpath(diry, "ovlp.2cluster.distance.subspace.yasashi.png")
    output = joinpath(diry, "ovlp.2cluster.distance.subspace.png")
    plot!(q, yoko, line1, type=:scatter, markershapes=:circle, label=labelname, legendfontsize = 15, tickfontsize = 15)
    plot!(r, yoko, line1, type=:scatter, markershapes=:circle, label=labelname, legendfontsize = 15, tickfontsize = 15, labelfontsize=15)
    savefig(r, output)
    for j in axes(darray,1)
        offset = 0.05 # オフセットの値を調整する
        annotate!(q, yoko[j], line1[j] + offset, text(string(round(line1[j], digits=3)), 8))
    end
    savefig(q,outputgentle)
    display(q)
end




#R<=10配位
Rtrn = 10
gatefull = fill(true, length(gate))
Rtrngate = gatefull
for i in 1:length(Rtrngate)
    if coord[i,1] > Rtrn
        Rtrngate[i] = false
    end
end 
subspace_plot("R$(Rtrn)trn", Rtrngate, 1.0e-15)
Hvec_trn
norvec
#θ=>π/4
θtrn = π/4
subname ="θπ4trn"
begin
    θtrngate = fill(true, length(gate))
    for i in 1:length(gatefull)
        if coord[i,3] < θtrn
            θtrngate[i] = false
        end
    end
    θtrngate = θtrngate .| gate6he .| clgate .| lingate
    trnt = 1.0e-15
    subspace_plot(subname, θtrngate, trnt)
end

#6he+2nのみを除く
#=
gateN6he = .!gate6he
subspace_plot("No6he", gateN6he, trn)
=#

#closure + 6he+2n + θ=90 配位
gate90s = gate90 .| clgate .| gate6he
subspace_plot("cl+6he+90", gate90s, 1.0e-15)

#suhara
suhara = gate90 .| clgate
subspace_plot("cl+90", suhara, trn)

#6he only
he690 = gate90 .| gate6he
subspace_plot("6he+90", he690, trn)

gateNols = gate6he .| clgate 
gateNols = .!gateNols
#ls抜く
subspace_plot("ls0", gateNols, trn)


#θ=90配位がメインでない、あるいは積分が必要であることを見るためにR:fixedでYのベクトル座標でオーバーラップをとる
for Rtrn in 3:12
    gatefull = fill(true, length(gate))
    Rtrngate = gatefull
    for i in eachindex(Rtrngate)
        if abs(coord[i,1] - Rtrn) > 0.1
            Rtrngate[i] = false
        end
    end 
    begin
        gate3body = clgate .| gate6he
        gate3body = .!gate3body
    end
    Rtrngate = Rtrngate .& gate3body
    #gateを使用してcoordから抜き出す
    αtrn = coord[Rtrngate,2]
    θtrn = coord[Rtrngate,3]

    Y = Rtrn .* sin.(αtrn)
    wfp1 = real.(wf[Rtrngate,1])
    wfp2 = real.(wf[Rtrngate,2])

    Yx = Y .* cos.(θtrn) ./ Rtrn
    Yy = Y .* sin.(θtrn) ./ Rtrn

    # 補間を行う新しいグリッドの生成
    lingrid_y = LinRange(minimum(Yy), maximum(Yy), 100)
    lingrid_x = LinRange(minimum(Yx), maximum(Yx), 100)
    grid_x, grid_y = np.meshgrid(lingrid_x, lingrid_y)

    # griddataを使用して補間し、範囲外の点は0で埋める
    zi = sc_interpolate.griddata((Yx, Yy), wfp1, (grid_x, grid_y), method="cubic", fill_value=0)

    output = joinpath(diry, "Rfixplot/R$(Rtrn)st1.png")
    p = contourf(lingrid_x, lingrid_y, zi, xlabel="Yx", ylabel="Yy", xrange=(0,1.0), yrange=(0,1.0), aspect_ratio=:equal, c=:thermal, title="R=$(Rtrn)")
    savefig(p, output)

    # griddataを使用して補間し、範囲外の点は0で埋める
    zi = sc_interpolate.griddata((Yx, Yy), wfp2, (grid_x, grid_y), method="cubic", fill_value=0)
    output = joinpath(diry, "Rfixplot/R$(Rtrn)st2.png")
    p = contourf(lingrid_x, lingrid_y, zi, xlabel="Yx", ylabel="Yy", xrange=(0,1.0), yrange=(0,1.0), aspect_ratio=:equal, c=:thermal, title="R=$(Rtrn)")
    savefig(p, output)
end
evol
evolls0 
norvec
if evol == evolls0
    println(true)
end
N = 490
trnum = N-tN +1
#ls力0の解析
#suharaにあるように、固有エネルギーのlsに対する変化
Etrn
Etrnls0
binary = "./$diry/resiroiro.jld2"
cls = GCM_keisu(evolls0, Hvec_trnls0, norvecls0)
begin
    Ewols = Vector{Vector{Float64}}()
    keisuls = Vector{Matrix{ComplexF64}}()
    push!(Ewols,Etrnls0)
    cls = GCM_keisu(evolls0, Hvec_trnls0, norvecls0)
    push!(keisuls, cls)
    HΔls = H - Hls0
    HΔls = HΔls ./ 1600
    gatefull = trues(N)
    #ls400刻み0から3200まで
    for ls in 200:200:3200
        Hls = Hls0 + HΔls .* ls
        Hlsvec, evolls, norvecls, Etrnls = calcu_matome(gatefull, Hls, ol, trn)
        cls = GCM_keisu(evolls,Hlsvec,norvecls)
        push!(Ewols, Etrnls)
        push!(keisuls, cls)
    end

    binary = "./$diry/resiroiro.jld2"
    @save binary Ewols keisuls
end
#Els[16]
#@load binary Els

Hls = zeros(Float64,Npick,16)
yokols = zeros(Float64, 16)
for st in 1:Npick
    for j in 1:16
        Hls[st,j] = Ewols[j][st]
        if j > 1
            yokols[j] = 400+(j-2)*200
        end
    end
end

for st in 1:Npick
    Elsp = Hls[st,:]
    if st==1 #出力の絵を１枚にする。
        global p = plot(xlabel="ls parameter", ylabel="Energy[MeV]")
        labelname = "g.s."
    elseif st==2 
        labelname = "2nd"
    else 
        labelname = "$(st)th"
    end
    plot!(p, [yokols], [Elsp], type=:scatter, markershapes=:circle, label=labelname)
end

display(p)
output = "./$(diry)/lsvary.png"
savefig(p,output)