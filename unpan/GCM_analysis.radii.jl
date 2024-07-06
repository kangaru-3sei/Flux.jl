#平均二乗半径を.jld形式で出力するだけ


#以下read_dataからのコピペ
using JLD2, LinearAlgebra
begin
    #識別子を読み込む
    #println("解析の識別子を打ち込んでね")



    #ファイル保存場所
    Rdiry = "./$(Rname)"
    isdir(Rdiry) || prinrln("フォルダ名の打ち間違いだよ！")

    #アウトプットのバイナリ。カレントディレクトリに出力
    Rbinary = "radii.$(Rname).jld2"
    Rbinary = joinpath(Rdiry, Rbinary)

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

    H_vec,overlap_vec, N = read_data(file)

    H_mat = vec_mat(H_vec)
    #H_mat = Hermite(H_mat)
    H_mat = Hermitian(H_mat)
    #H_mat = Symmetric(H_mat)

    return H_mat, N
end

#処理
begin
    file = joinpath(Rdiry, "temps5.$(Rname).out")
    radius, N= eigen_hamiltonian(file)
#バイナリにつっこむ
    @save Rbinary radius
    println("行列要素のサイズは$N")
    println("ほらよ解析しな")
    println("radii.$(Rname).jld2にr.m.s. radii の行列を入れとる")
    #println("subspaceのovlpを見たいがためにnorm kernelの打ち切りtrn、規格化されたノルムカーネルの固有値evol固有ベクトルnorvec、truncateされて残った数tN、proHの固有ベクトルHvec_trnを記録しとる")



    #=
    for i in 0:num_state
        overlap_figure(result, i, diry)
    end
    =#
end