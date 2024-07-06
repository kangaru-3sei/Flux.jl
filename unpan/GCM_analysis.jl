#Ring-SchuckにあるGCM計算の解析手法を用いる。
#ノルムカーネルの固有値方程式を解き、一定以下の固有値のtruncationを行う。
#射影したヒルベルト空間における固有値方程式を解く。

#各配位のd1, d2,θ の値、指定した状態とのオーバーラップの値をresultで出力


#以下read_dataからのコピペ
using  LinearAlgebra, JLD2

begin
    #閾値
    trn = 3.0e-16

    #オーバーラップをとる状態の数
    Npick = 8
end

begin
    #識別子を読み込む
    println("解析の識別子を打ち込んでね")
    name = readline()


    #ファイル保存場所
    diry = "./$(name)"
    isdir(diry) || mkdir(diry)

    #このプログラムのinput
    filename = "temps5.$(name).out"
    file = joinpath(diry, filename)
    isfile(file) || cp(filename, file)

    inputname = "temps5.$(name).inp"
    input = joinpath(diry, inputname)
    isfile(input) || cp(inputname, input)

    inpcoor = joinpath(diry, "input.coor.out")

    #アウトプットファイル
    output = "./$(name)/ovlp.$(name).out"

    #アウトプットのバイナリ。カレントディレクトリに出力
    binary = "Eovl.$(name).jld2"
    binary = joinpath(diry, binary)

    #inpbinary = "./$(name)/coor$(name).jld2"

    #hedderread = "input.coor.out"
    #hedderread = joinpath(diry, hedderread)

end



function eigen_hamiltonian(file)
#ファイルから数字を読み取り、行列要素と重なり積分をベクトルとして出力。
    function read_data(file)
        nind = 0  # 行列の独立な成分の数。jgiveの個数を数える。
        nm = 0  #出力された各値をベクトルに収納するさいに、ベクトルの成分番号を指定する。
        N = 0 #行列のサイズをnmから出す。下に記述。
        gate = false  # ゲートの開閉によりjgiveの個数を数えたりjgive以下の二個のベクトルを別別に読み込んだりさせる（詳細は下）。
        gate2 = false   #ハミルトニアンに入れるか、重なり積分に入れるかを判別する

        open(file, "r") do file
            for line in eachline(file)  #最初に独立な行列要素の数を勘定
                if occursin("jgive", line)
                    nind += 1
                end
            end
        end

        #格納するための箱を定義
        H_vec = zeros(Complex{Float64}, nind)  # ハミルトニアンの行列要素
        #H_vec = zeros(Float64, nind)
        overlap_vec = zeros(Complex{Float64}, nind)  # 波動関数の重なり積分
        #overlap_vec = zeros(Float64, nind)

        #中身を書き出す
        open(file, "r") do file
            for line in eachline(file)
                if occursin("jgive", line)
                    gate = true
                    nm += 1  # 独立な行列要素の数
                elseif gate == true
                    if occursin("(", line)
                        num_str = split(replace(line, r"[)(]" => ""), ",")
                        num = complex(parse(Float64, num_str[1]), parse(Float64, num_str[2]))
                        #num = parse(Float64, num_str[1])

                        if gate2 == false  # 数字の代入。H_matが空かどうかで判断。
                            H_vec[nm] = num
                            gate2 = true    #次の数字は重なり積分のほうに入れる。
                        else
                            overlap_vec[nm] = num
                            gate = false    # 重なり積分に収納したら、ゲートを閉じる。
                            gate2 = false   #次の数字は行列要素のほうに入れることにする。
                        end
                    end
                end
            end
        end

        #行列の独立な成分の数nmから行列サイズNを出す。
        while nm > 0
            N += 1
            nm = nm - N
        end
        println("アウトプットファイルから読み取った行列サイズ：",N)
        return H_vec, overlap_vec, N
    end


    # ベクトルを片半分行列に直す関数
    function vec_mat(vec)
        mat = zeros(Complex{Float64}, N, N)
        #mat = zeros(Float64, N, N)
        for i in 1:N
            for j in i:N
                mat[i, j] = vec[1]
                deleteat!(vec, 1)
            end
        end
        return mat
    end


    # エルミート行列に変換する関数
    function Hermite(mat)
        for i in 1:N
            for j in 1:N
                if i > j
                    mat[i, j] = conj(mat[j, i])
                elseif i == j
                    mat[i,j] = real(mat[i,j])
                end
            end
        end
        return mat
    end

    H_vec, ol_vec, N = read_data(file)

    H_mat = vec_mat(H_vec)
    #H_mat = Hermite(H_mat)
    H_mat = Hermitian(H_mat)
    #H_mat = Symmetric(H_mat)

    ol_mat = vec_mat(ol_vec)
    #ol_mat = Hermite(ol_mat)
    ol_mat = Hermitian(ol_mat)
    #ol_mat = Symmetric(ol_mat)

    #一般化固有値問題を解く
    #ε,ψ = eigen(H_mat,ol_mat)

    evol, ecol = eigen(ol_mat)

    return H_mat, ol_mat, N, evol, ecol
end



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


#trnol = integ_eivec(norvec, trnum, N)

#ノルムカーネルが基本行列になっているか確認
function truned_OL(evol, norvec, trnum, ol, N)
    tN = N - trnum + 1      #projectionを行った後の行列のサイズ
    proOL = zeros(Complex{Float64}, tN, tN)
    for k in 1:tN       #行列要素の成分番号
        re_k = k + trnum - 1    #元の固有ベクトル行列に対応する横の番号
        for l in k:tN       #行列要素の成分番号
            proOL_elem = 0
            re_l = l + trnum - 1
            #生成座標による積分
            for i in 1:N    
                for j in 1:N
                    proOL_elem += conj(norvec[i,re_k])/sqrt(evol[re_k])*norvec[j, re_l]/sqrt(evol[re_l])*ol[i,j]
                end
            end
            proOL[k,l] = proOL_elem
        end
    end
    proOL = Hermitian(proOL)
    return tN, proOL
end

#=
tN, trnOL = truned_OL(evol, norvec, trnum, ol, N)
trnOL
=#


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
end


#自然基底を用いたoverlapの計算
function exol(basenum, Hvec_trn, ol, evol, norvec, trnum, tN, N, state_num)
    result = 0.0+0.0im
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

using Printf
using Plots
#インプットからd1, d2, θを読み取りoverlap期待値とともに記録
#さらに各配位の空間位置を記録する
function find_degenerate(input, inpcoor, output, H, ol, evol, norvec, trnum, Hvec_trn, tN, N, num_state)
    ν = 0.5 #適当に定義

    #現在の行数
    current_line = 1

    #核子数を読み出すところの行数
    line_nuc = 1

    #νを読み出すところの行数
    line_ν = 4

    #各配位のエネルギー期待値とオーバーラップ。後でまとめてバイナリに出力
    Earray = zeros(Float64, N)
    ovlarray = zeros(Complex{Float64}, N, num_state)


    open(input, "r") do input
        for line in eachline(input)
            #核子数を読み出す
            if current_line == line_nuc
                #空白で分割して数値の配列に変換
                numbers = split(line) |> x -> parse.(Int16, x)
                
                pnum = numbers[1] + numbers[2]
                nnum = numbers[3] + numbers[4]
                global num_nuc = pnum + nnum 
            end
            #波束幅を読こむ
            if current_line == line_ν
                #空白で分割して数値の配列に変換
                numb = split(line) |> y -> parse.(Float64, y)
                # 4行目の最初の数値を取得
                ν = numb[1]
                break
            end
            current_line += 1
        end
    end

    #パラメータ部分をとばす
    line_num = 7

    #current_lineを戻す
    current_line = 1

    gate = 1        #x,y,z座標を順番に読み取る
    coor = zeros(Complex{Float64}, 3, num_nuc)
    #g = zeros(Complex{Float64}, 3)
    particle_num = 1   #粒子数を判別するためのゲート
    
    #読み取った配位の数を確認
    n = 1       


    open(output, "w") do output
        #アウトプットのヘッダー
        #input.coorから読み取る
        open(inpcoor, "r") do inpcoor
            counter = 1
            for line in eachline(inpcoor)
                if counter == n
                    write(output, line)
                    break
                end
                counter += 1
            end
        end

        print(output, ",E")
        for i in 1:num_state
            print(output, ",ol", i)
        end
        println(output)
    end

    open(input, "r") do input
        for line in eachline(input)
            if current_line == line_num
                if gate == 1
                    num_str = split(replace(line, r"[)(]" => ""), ",")
                    global coor[1, particle_num] = parse(Float64, num_str[1]) + parse(Float64, num_str[2])*im
                    gate += 1     #次に読み取る数字はy座標
                    line_num += 1   
                #y座標を読み取る
                elseif gate == 2
                    num_str = split(replace(line, r"[)(]" => ""), ",")
                    global coor[2, particle_num] = parse(Float64, num_str[1]) + parse(Float64, num_str[2])*im
                    gate += 1
                    line_num += 1
                #z座標を読み取る
                elseif gate == 3
                    num_str = split(replace(line, r"[)(]" => ""), ",")
                    global coor[3, particle_num] = parse(Float64, num_str[1]) + parse(Float64, num_str[2])*im
                    gate = 1
                    line_num += 1
                    particle_num += 1
                end

                if particle_num == num_nuc + 1
                    open(output, "a") do output
                        #幾何学的配置を算出
                        #内積から角度を計算
                        
                            #三体の場合のhyperradius
                            
                        #ベクトルR1,　R2とその間のなす角を計算
                        #=
                        R1 = [x1-x2, z1-z2] * √ν
                        CM2n = [x1+x2, z1+z2] * 0.50 
                        R2 = [CM2n[1]-xα, CM2n[2]-zα] * √ν
                        R1nor = norm(R1)
                        R2nor = norm(R2)
                        R = sqrt((R1nor)^2 + 2*(R2nor)^2)
                        cos_theta = dot(R1, R2)/ (norm(R1)*norm(R2))
                        #acos関数に渡す前に補正
                        clamp(x) = max(min(x, 1.0), -1.0)
                        θ = acos(clamp(cos_theta)) * 180/π
                        =#


                        #各配位のエネルギー期待値
                        Earray[n] = H[n,n]/ol[n,n]
                        E = Earray[n]

                        # input.coor.outの同じ行をoutputに持ってくる
                        open(inpcoor, "r") do inpcoor
                            counter = 1     
                            for line in eachline(inpcoor)
                                if counter == n+1   #hedderは除く
                                    write(output, line)
                                    break
                                end
                                counter += 1
                            end
                        end

                        # フォーマットされた文字列を作成してファイルに書き込む
                        # 各列の幅を指定
                        header = @sprintf(",%.2f", E)
                        write(output, header)

                        # 指定した状態の数だけオーバーラップを記録
                        for j in 1:num_state
                            expedol = exol(n, Hvec_trn, ol, evol, norvec, trnum, tN, N, j)
                            ovlarray[n,j] = expedol
                            detail = @sprintf(",%.2f", real(expedol))
                            write(output, detail)
                        end

                        #重心を記録
                        #=
                        g = sum(coor, dims=2)
                        for z in g
                            @printf(output, "%.2f+%.2fim ", real(z), imag(z))
                        end
                        =#
                        #kakuninファイルを改行
                        println(output)
                    
                    end

                    #配位の順番数を更新
                    n += 1
                    particle_num = 1
                    
                    
                    #スピン部分を飛ばす
                    line_num += num_nuc*2 + 2 
                end
            end
            current_line += 1
        end
    end
    return Earray, ovlarray
end



#処理
begin
    H, ol, N, evol, evcol = eigen_hamiltonian(file)
    norvec, trnum = normed_eivec(evol, evcol, trn, N)
    proH, Etrn, Hvec_trn, tN = truned_H(evol, norvec, trnum, H, N)
    println("固有エネルギーを記録の個数分だけ")
    tN = N - trnum + 1
    println("射影空間における基底の数：", tN)
    num_state = Npick > tN ? tN : Npick 
    for i in 1:num_state
        println(Etrn[i])
    end
    E, wf = find_degenerate(input, inpcoor, output, H, ol, evol, norvec, trnum, Hvec_trn, tN, N, num_state)

    #norvec,evolがなぜか使用しない固有ベクトルが残ったままだったので切る。
    norvec = norvec[:,trnum:end]
    evol = evol[trnum:end]

        #input.coor.outが複製されたので消す
        #rm(hedderread)
    
#バイナリにつっこむ
    @save binary proH Etrn wf H ol trn evol norvec tN Hvec_trn
    println("ほらよ解析しな")
    println("Eovl.$(name).jld2に無加工H、無加工ノルムカーネルol、自然基底射影したproH、各配位のエネルギーEtrn、オーバーラップwfを記録しとる。")
    println("subspaceのovlpを見たいがためにnorm kernelの打ち切りtrn、規格化されたノルムカーネルの固有値evol固有ベクトルnorvec、truncateされて残った数tN、proHの固有ベクトルHvec_trnを記録しとる")



    #=
    for i in 0:num_state
        overlap_figure(result, i, diry)
    end
    =#
end