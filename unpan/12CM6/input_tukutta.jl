#12Cの場合の、hyperradius(R,α)とR₁、R₂のなす角θによるGCM配位の作成。
#Rさまざま、full配位。500以下に抑える
#角度方向の基底数はRの大きさによって変わる
#interactionは比較のためsuharaを採用

begin
    #あと、フラグ管理でshell配位を入れるかを選べる
    Np = 6
    Nn = 6
    M = 0.60    #Majorana term
    B = 0.0 #Bartlett term
    H = B   #Heisenberg term
    v1 = 1600  #ls項の係数其一
    v2 = -v1    #ls項の係数その2
    ν = 0.235    #ガウス幅パラメータ
    Jmax = 0
    flag_shell = true  #shell配位を入れるかを選ぶ

    name = "12CM6"

    Rmax = 12

    αmax = π/2
    αmin = 0
    thr1 = π/2/19
    function dα(R,thr)
        d = π/2/R
        if d>thr
            return d
        else 
            return thr
        end
    end
    
    θmax = π/2
    θmin = 0
    thr = π/2/10
    function dθ(R,α, thr)
        d = π/2/(R*sin(α))
        if d>thr
            return d
        else 
            return thr
        end
    end

    d = 0.0    #クラスター間距離がこれより小さくなったら処理を行わない。
    
    #name = readline()
    
end

include("./makeconfig.jl")
include("./AQCMconfig.jl")
include("./3bodyconfig.jl")

function posit12C_gate(R, α, θ, d, ν)
    #原点を重心にもってきた３体配位の作成
    R1 = R*cos(α)/√2
    R2 = R*sin(α)*sqrt(3/8)
    x1 = 2/3*R2*sin(θ)*√ν
    x = -1/3*R2*sin(θ)*√ν
    z1 = 1/3*R2*cos(θ)*√ν
    z2 = (0.50*R1-1/6*R2*cos(θ))*√ν
    z3 = (-0.50*R1-1/6*R2*cos(θ))*√ν

    #２クラスター間の距離が近すぎたら計算をしない
    D12 = sqrt((z1-z2)^2+(x1-x)^2)
    D23 = R*cos(α)
    
    if D12 > d && D23 > d
        return true
    else 
        return false
    end
end

#αの刻み幅を決める
#=
begin
    R = 14
    mx=0
    n=19
    #for n in 10:20
        for i in 0:n-1
            dx8 = R*(cos((i+1)/n/2*π)-cos(i/n/2*π))
            dx8 = abs(dx8)
            dy8 = R/√2*(sin((i+1)/n/2*π)-sin(i/n/2*π))
            dy8 = abs(dy8)
            dx12 = R/√2*(cos((i+1)/n/2*π)-cos(i/n/2*π))
            dx12 = abs(dx12)
            dy12 = R*√(3/8)*(sin((i+1)/n/2*π)-sin(i/n/2*π))
            dy12 = abs(dy12)
            mx = max(mx,dx8,dy8,dx12,dy12)
        end
        #if mx<1.3
        #    break
        #end
    #end
    #dr = R/√2*π/2/m
    #mx = max(mx,dr)
    println(mx)
end
=#


#=
for i in 1:14
    R = i
    Nbase = Nfunc(R, flag_shell, flag_6he)
    println(Nbase)
end
=#

function HSP12C(n,sakusei)
    # 初めにフラグが立っているときシェル的な配位を入れる。
    if flag_shell == true
        n = C12closure(ν,n)   
    end

    for R in 1:Rmax
        #R=1の配位が直線配位しかなく３体配位とのovlpが0になったので変更
        if R > 1
            #Y=0を先に加えておく
            α = 0
            θ = 0
            n = α3body(R, α, θ, d, ν, n)

        else
            #Y=0を先に加えておく
            α = π/4
            θ = π/2

            n = α3body(R, α, θ, d, ν, n)

        end


        #Yが１以上の配位を加える
        α = dα(R,thr1)  #α=0は含まないのであらかじめ足しておく
        while α<π/2-0.001
            θ = 0
            while θ <= π/2 - dθ(R, α, thr)/10   #θ=π/2が近いところまで。必ずθ=π/2を含めるため、近くは除く

                if posit12C_gate(R, α, θ, d, ν) == true
                    n = α3body(R, α, θ, d, ν, n)  
                    θ += dθ(R,α,thr)
                end          

            end

            θ = θ-dθ(R,α,thr)    #判定のため引き戻し
            #θ=π/2が含まれないときは含める
            if θ != π/2
                θ = π/2
                n = α3body(R, α, θ, d, ν, n)                    
                
            end
            α += dα(R,thr1)  
        end
    end

    sakuseigo_atoshori(hedder, kitei_kakunin, sakusei, Np, Nn, Jmax, M, B, H, v1, v2, ν, n)
end
HSP12C(n,PROGRAM_FILE)