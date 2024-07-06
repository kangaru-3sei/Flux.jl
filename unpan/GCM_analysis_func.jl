#GCM_analysis.jlの関数だけ
#subspaceにおけるbareなH, bareなol, trnを代入



#固有ベクトルの規格化
function normed_eivec(evol, evcol, trn, N)
    trnum = 1
    for i in 1:N-1
        if evol[i] <= trn && evol[i+1] >= trn
            trnum = i+1   #ここを含むこれより先のみ取る、という閾値の固有値の番号。
        end
    end
    #固有ベクトルの規格化
    norvec  = zeros(Complex{Float64}, N, N)
    for i in trnum:N
        vecsum = 0
        for j in 1:N
            vecsum += abs(evcol[j,i])^2
        end
        vec_norm = sqrt(vecsum)
        for j in 1:N
            norvec[j,i] = evcol[j,i]/vec_norm
        end
    end
    return norvec, trnum
end

#固有ベクトルが直交しているか確認
function integ_eivec(norvec, trnum, N)
    tN = N - trnum + 1 
    trnol = zeros(Complex{Float64}, tN, tN)
    for k in 1:tN
        re_k = k + trnum - 1 
        for l in k:tN
            trnol_elm = 0
            re_l = l + trnum - 1
            for i in 1:N
                trnol_elm += conj(norvec[i, re_k])*norvec[i,re_l]
            end
            trnol[k,l] = trnol_elm
        end
    end 
    return trnol
end    


#射影空間におけるハミルトニアン
function truned_H(evol, norvec, trnum, H, N)
    tN = N - trnum + 1      #projectionを行った後の行列のサイズ
    proH = zeros(Complex{Float64}, tN, tN)
    for k in 1:tN       #行列要素の成分番号
        re_k = k + trnum - 1    #元の固有ベクトル行列に対応する横の番号
        for l in k:tN       #行列要素の成分番号
            proH_elem = 0
            re_l = l + trnum - 1
            #生成座標による積分
            for i in 1:N    
                for j in 1:N
                    proH_elem += conj(norvec[i,re_k])/sqrt(evol[re_k])*norvec[j, re_l]/sqrt(evol[re_l])*H[i,j]
                end
            end
            proH[k,l] = proH_elem
        end
    end
    proH = Hermitian(proH)
    Etrn, Hvec_trn = eigen(proH)
    return proH, Etrn, Hvec_trn, tN
end

#状態ベクトルの規格化
function norm_state(Hvec_trn, tN, state_num)
    result = 0
    for i in 1:tN
        result += abs2(Hvec_trn[i,state_num])
    end
    result = √result
    result
end


#自然基底を用いたoverlapの計算
function exol(basenum, Hvec_trn, ol, evol, norvec, trnum, tN, N, state_num)
    result = 0
    for k in 1:tN   #自然基底の個数
        re_k = k + trnum - 1    #元の固有ベクトル行列に対応する横の番号
        sum = 0 
        for i in 1:N    #生成座標による総和
            sum += norvec[i,re_k]*ol[basenum,i]
        end
        result += Hvec_trn[k,state_num]/√evol[re_k]*sum
    end
    result = result/√ol[basenum, basenum]
end

#gateで指定した成分を除去
#subspaceとのovlpを計算するため
function matremove(mat, mask)
    #matとgateのサイズが違う場合、警告を鳴らす
    if length(mask) != size(mat,1)
        println("gateと行列のサイズが違いますわ...")
        return 0
    end
    # gateのfalseの位置を見つける
    line_remove = findall(.!mask)

    # 指定した列を除去
    mat = mat[:, setdiff(1:end, line_remove)]
    # 指定した行を除去
    mat = mat[setdiff(1:end, line_remove), :]

    #hermiteにする
    mat = Hermitian(mat)
    return mat
end
#gateで指定したsubspaceとのovlpを計算するため
function calcu_matome(gate, bareH, bareol, trn)
    subH = matremove(bareH, gate)
    subol = matremove(bareol, gate)

    evol, evcol = eigen(subol)
    N = size(subol,1)
    norvec, trnum = normed_eivec(evol, evcol, trn, N)
    proH, Etrn, Hvec_trn, tN = truned_H(evol, norvec, trnum, subH, N)
    tN = N - trnum + 1
    #norvecがなぜか切り捨て固有ベクトルを残した形になっていたので変更
    norvec = norvec[:,trnum:end]
    #ovlの計算に使用するevolも変更
    evol = evol[trnum:end]
    return Hvec_trn, evol, norvec, Etrn
end

function wfsubspace(bareH, bareol, trn)
    evol, evcol = eigen(bareol)
    N = size(bareol,1)
    norvec, trnum = normed_eivec(evol, evcol, trn, N)
    proH, Etrn, Hvec_trn, tN = truned_H(evol, norvec, trnum, bareH, N)
    tN = N - trnum + 1
    println("元々の行列要素の個数：", N)
    println("射影空間における基底の数：", tN)
    num_state = Npick > tN ? tN : Npick 
    println("固有エネルギーを記録の個数分だけ")
    for i in 1:num_state
        println(Etrn[i])
    end
    #ovl
    ovlarray = zeros(Complex{Float64}, N, num_state)
    for n in 1:N
        for j in 1:num_state
            expedol = exol(n, Hvec_trn, bareol, evol, norvec, trnum, tN, N, j)
            ovlarray[n,j] = expedol
        end
    end
    return ovlarray
end

#subspaceとのovlpを撮るときに使用するHvec_trn, evol, norvecを返す。
function Hvec_trn_sub(bareH, bareol, trn)
    evol, evcol = eigen(bareol)
    N = size(bareol,1)
    norvec, trnum = normed_eivec(evol, evcol, trn, N)
    proH, Etrn, Hvec_trn, tN = truned_H(evol, norvec, trnum, bareH, N)
    return Hvec_trn, evol, norvec
end

#HvecをbareHのサイズに戻す
function Hvec_kakucho(Hvec_sub, n::Int)
    N = size(Hvec_sub,1)
    A = zeros(ComplexF64, N, n-N)
    res = hcat(A,Hvec_sub)
    A = zeros(ComplexF64, n-N, n)
    res = vcat(A,res)
    return res
end

function matremove(mat, mask)
    #matとgateのサイズが違う場合、警告を鳴らす
    if length(mask) != size(mat,1)
        println("gateと行列のサイズが違いますわ...")
        return 0
    end
    # gateのfalseの位置を見つける
    line_remove = findall(.!mask)

    # 指定した列を除去
    mat = mat[:, setdiff(1:end, line_remove)]
    # 指定した行を除去
    mat = mat[setdiff(1:end, line_remove), :]

    #hermiteにする
    mat = Hermitian(mat)
    return mat
end

function ovl_sub_ochiai(H,ol,mask, gate, Npick, name, diry, coord, trn)
    H = matremove(H, mask)
    ol = matremove(ol, mask)
    wf = wfsubspace(H, ol, trn)

    #w.f.をテクストファイルで出力
        coordp = coord[mask,:]  #まずsubspaceの座標を抜き出す
        X = coordp[:,1] .* cos.(coordp[:,2])
        Y = coordp[:,1] .* sin.(coordp[:,2])
        θ = coordp[:,3] .* 180 ./ π
        
        #subspaceから3bodyの配位を選ぶ
        gatecoord = gate[mask]

    output = joinpath(diry, "ovlp$name.txt")
    open(output, "w") do file
        write(file, "n 3体 X Y θ")
        for i in 1:Npick
            write(file, " st$i")
        end
        println(file)
        for base in 1:size(H,1)
            write(file, "$base ")
            #3bodh gate
            if gatecoord[base] == true
                write(file, "T ")
            else 
                write(file, "F ")
            end

            #座標
            str = @sprintf("%.2f %.2f %.1f ", X[base], Y[base], θ[base])
            write(file, str)

            #wave function
            for st in 1:Npick
                str = @sprintf("%.3f ", real(wf[base, st]))
                write(file, str)
            end
            println(file)
        end
    end

    #plotに使用する条件を満たす点を渡す
    gate3b = gate .& mask   
    masksubθ = coord[:,3] .== π/2 
    masksubY0 = coord[:,2] .== 0
    masksub = masksubY0 .| masksubθ     
    gateplot = gate3b .& masksub    #全490配位の中で、subspaceかつd1=d2のプロット条件を満たすもの
    coordp = coord[gateplot,:]    #次にsubspaceの中から空間座標としてプロットする者を抜き出す。

    gatewfp = gateplot[mask]
    wf = wf[gatewfp,:]
    return wf, coordp
end

ovl_sub(mask, name, Npick, diry, trn) = ovl_sub_ochiai(H, ol, mask, gate, Npick, name, diry, coord, trn)

#subspaceとのオーバーラップをθ=π/2の配位でとる
function XYedplot(diry, name, mask, Npick, trn)
    println("********$name********")
    wfeqd, coeqd = ovl_sub(mask, name, Npick, diry, trn)
    #Y,X/2の座標でプロット
    Yeqd = coeqd[:,1] .* sin.(coeqd[:,2])
    Xeqd = coeqd[:,1] .* cos.(coeqd[:,2]) 
    #いつもの
    for i in 1:Npick    # ith state
        wfp = abs2.(wfeqd[:,i])    
        
        # 補間を行う新しいグリッドの生成
        lingrid_y = LinRange(minimum(Yeqd), maximum(Yeqd), 100)
        lingrid_x = LinRange(minimum(Xeqd), maximum(Xeqd), 100)
        grid_x, grid_y = np.meshgrid(lingrid_x, lingrid_y)
        
        # griddataを使用して補間し、範囲外の点は0で埋める
        zi = sc_interpolate.griddata((Xeqd, Yeqd), wfp, (grid_x, grid_y), method="cubic", fill_value=0)

        output = joinpath(diry, "wf$(name).st$i.png")
        p = contourf(lingrid_x, lingrid_y, zi, xlabel="X", ylabel="Y", xrange=(0,12), yrange=(0,12), clims=(-0.01, 0.45), aspect_ratio=:equal, c=:thermal, title="ovlp. of $i th state")
        savefig(p, output)
    end
end

subspace_plot(name, mask, trn) = XYedplot(diry, name, mask, 5, trn)

function GCM_keisu(evol,Hvec_trn,norvec)
    tN = length(evol)
    N = size(norvec,1)
    #gcm配位につく係数を求める: ith eigstate |Ψ>=Σ_j c[j,i] ϕ(j)
    c = Matrix{Complex{Float64}}(undef, N, tN)
    for j in 1:N
        for i in 1:tN
            c[j,i] = sum(Hvec_trn[k,i]/√evol[k]*norvec[j,k] for k in 1:tN)
        end
    end
    return c
end

#subspaceとのオーバーラップ。subspaceでのnatural basisとj th state in full spaceとのオーバーラップの値を出力 
#inputはgate:subspace指定 のみ。
function calcu_subovl(gate, bareH=H, bareol=ol, eigvals=evol, threshold=trn, proHvec=Hvec_trn, provec=norvec)
    subHpro_vec, subolevals, subnorvec, Esubtrn = calcu_matome(gate, bareH, bareol, threshold)
    
    C = GCM_keisu(eigvals,proHvec,provec)
    c = GCM_keisu(subolevals, subHpro_vec, subnorvec)
    
    output = Vector{Float64}(undef, 8)  #natural basisとのovlを保存。多くても面倒なので第8状態まで
    
    function calcu_sqrdovlp(i, F=C, f=c, N=bareol, subgate=gate)
        line_remain = findall(subgate)
        conjF = conj.(F)
        res = sum(abs2(sum(conjF[a,i]*f[α,j]*N[a,line_remain[α]] for a in axes(F,1) for α in axes(f,1))) for j in axes(c,2))
        return res
    end

    for i in axes(output,1)
        output[i] = calcu_sqrdovlp(i)
    end

    return output
end

#ある物理量の期待値を求める。状態は固有所帯iの fullのみ
#radiusは求めたい物理量をgcm配位で挟んだもの
function expecval_full(radius::Hermitian,i,F=c)
    N = size(F,1)
    #radiusの期待値
    res = sqrt(real(sum(conj(F[k,i])*F[l,i]*radius[k,l] for k in 1:N, l in 1:N)))
    
    return res
end