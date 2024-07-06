#Ring-SchuckにあるGCM計算の解析手法を用いる。
#ノルムカーネルの固有値方程式を解き、一定以下の固有値のtruncationを行う。
#射影したヒルベルト空間における固有値方程式を解く。

#shell配位を除いたものを出力
#=これに伴う変更点
    1.N→N-5: line 148
    2.H_mat = H_mat[6:end, 6:end]: line 140(ovlも)
    2.line_num = 217: line 314
=#

#coor$(name).jld2から座標を読み取っているけど配列の名前に依存しているので注意(line 299,)

begin
    using LinearAlgebra
    using Glob
    using JLD2
end

begin
    #識別子を読み込む
    println("解析の識別子を打ち込んでね")
    name = readline()
end

begin
    #アウトプットのファイル数を読み取る
    function countfile(name)
        # 指定ディレクトリのファイルをリストアップ
        files = glob("./$(name)/temps5.$(name)*.inp", ".")

        # ファイルの数をカウント
        count = length(files)
        return count
    end

    n = countfile(name)
    println("ファイル数は$(n)ね")
end


begin
    #閾値
    global trn = 1.0e-8

    #オーバーラップをとる状態の数
    global Npick = 1
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

    #shell配位を除く
    #H_mat = H_mat[6:end, 6:end]

    ol_mat = vec_mat(ol_vec)
    #ol_mat = Hermite(ol_mat)
    ol_mat = Hermitian(ol_mat)
    
    #shell配位を除く
    #ol_mat = ol_mat[6:end, 6:end]

    evol, ecol = eigen(ol_mat)
    #N = N-5

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
    result = abs2(result)
end

using Printf
using Plots
#さらに各配位の空間位置を記録する
function find_degenerate(i, input, inpbinary, hedderread, output, H, ol, evol, norvec, trnum, Hvec_trn, tN, N, num_state)
    ν = 0.5 #適当に定義

    #現在の行数
    current_line = 1

    #核子数を読み出すところの行数
    line_nuc = 1

    #νを読み出すところの行数
    line_ν = 4

    #各配位のエネルギー期待値とオーバーラップ。後でまとめてバイナリに出力
    Earray = zeros(Float64, N)
    ovlarray = zeros(Float64, N, num_state)

    #input.coor.outからヘッダーだけ読み取る
    function Read1(hedderread)
        open(hedderread, "r") do file
            hedder = readline(file)
            return hedder
        end
    end

    #hedder = Read1(hedderread)

    #input座標データを読み取る
    @load inpbinary HR gate2

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

    #第5配位から開始
    #line_num = 217
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
        print(output, hedder)
        for i in 1:num_state
            print(output, " ol", i)
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


                        # フォーマットされた文字列を作成してファイルに書き込む
                        # 各列の幅を指定
                        judge = gate2[i][n]==true ? "T" : "F"
                        coor1 = HR[i][n,1]
                        coor2 = HR[i][n,2]
                        header = @sprintf("%d %s %.3f %.3f %.3f", n, judge, coor1, coor2, E)
                        write(output, header)

                        # 指定した状態の数だけオーバーラップを記録
                        for j in 1:num_state
                            expedol = exol(n, Hvec_trn, ol, evol, norvec, trnum, tN, N, j)
                            ovlarray[n,j] = expedol
                            detail = @sprintf(" %.3f ", expedol)
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




begin
    dir = "$(name)"

        #アウトプットのバイナリ。カレントディレクトリに出力
        binary = "Eovl.$(name).jld2"
        binary = joinpath(dir, binary)

        #↑に格納するデータ。truncationした空間におけるエネルギー固有値とオーバーラップ
        Eout = Vector{Vector{Float64}}()
        ovlout = Vector{Matrix{Float64}}()
        proHout = Vector{Matrix{Complex{Float64}}}()

        #↑に格納するデータ2。元のHとノルムカーネル
        bareH = Vector{Matrix{Complex{Float64}}}()
        bareovl = Vector{Matrix{Complex{Float64}}}()

    for i in 1:n
        #ファイル保存場所
        diry = "./$(name)/$(i)"
        isdir(diry) || mkdir(diry)

        #このプログラムのinput
        filename = "temps5.$(name)$(i).out"
        file = joinpath(diry, filename)
        isfile(file) || cp(filename, file)

        inputname = "temps5.$(name)$(i).inp"
        input = joinpath(diry, inputname)
        isfile(input) || cp(inputname, input)

        inpbinary = "./$(name)/coor$(name).jld2"

        hedderread = "input.coor.out"
        hedderread = joinpath(diry, hedderread)

        #アウトプットファイル
        result = "ovlp.$(name)$(i).out"
        result = joinpath(diry, result)


        #処理
        H, ol, N, evol, evcol = eigen_hamiltonian(file)
        norvec, trnum = normed_eivec(evol, evcol, trn, N)
        proH, Etrn, Hvec_trn, tN = truned_H(evol, norvec, trnum, H, N)
        println("固有エネルギーを記録の個数分だけ")
        tN = N - trnum + 1

        for i in 1:Npick
            println(Etrn[i])
        end
        println("射影空間における基底の数：", tN)
        num_state = Npick > tN ? tN : Npick 
        #E, ovl = find_degenerate(i, input, inpbinary, hedderread, result, H, ol, evol, norvec, trnum, Hvec_trn, tN, N, num_state)

        #箱に突っ込む
        push!(Eout, Etrn)
        #push!(ovlout, ovl)
        push!(proHout, proH)
        push!(bareH, H)
        push!(bareovl, ol)

        #input.coor.outが複製されたので消す
        #rm(hedderread)
    end
#バイナリにつっこむ
@save binary Eout ovlout proHout bareH bareovl
println("ほらよ解析しな")
println("Eovl.$(name).jld2に無加工Hの行列要素bareH、無加工ノルムカーネルbareovl、自然基底でprojectionしたハミルトニアンproHout, 各配位のエネルギーEout、オーバーラップ期待値ovlを記録しとる。")
end


#=
for i in 0:num_state
    overlap_figure(result, i, diry)
end
=#