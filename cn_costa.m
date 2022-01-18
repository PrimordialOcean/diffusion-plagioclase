# Crank-Nicolson methodで実装
# 関数ファイルとして読み込まれるのを回避するためにダミーの値を入力
1;

# 関数定義
# CSVファイルを読み込む
function [distance, xan, wtmg] = get_data(filename = 'initial_value.csv')
    M = csvread(filename, 1, 0); # 先頭行をスキップ
    distance = M(:,1);
    xan = M(:,2);
    wtmg = M(:,3);
end

# 拡散係数を計算
# van Orman et al. (2014) を使用
function d = diff_coef(xan, tempc, dx, dt) # デフォルトは900℃
    tempk = tempc + 273.15;
    d = exp( - 6.06 - 7.96 * xan - 287000 / ( 8.31 * tempk) );
end

# 境界条件の部分を削除
function mat = cut_mat(mat)
    vecsize = size(mat);
    if vecsize(2) == 1
        mat(1) = [];
        mat(end) = [];
    else
        mat(1,:) = [];
        mat(end,:) = [];
        mat(:,1) = [];
        mat(:,end) = [];
    end
end

# 左辺側の係数行列を作成
function [matA, f_term1, f_term3] = calc_matA(nx, dt, dx, d, xan, tempc)
    # パラメータを指定
    aconst = -26100; # Bindeman et al. (1998)
    R = 8.314; # gas constant
    tempk = tempc + 273.15;
    # 係数を計算
    f_term1 = zeros(nx);
    f_term2 = zeros(nx);
    f_term3 = zeros(nx);
    for i = 2:nx-1
        f_term1(i) = ( d(i+1) - d(i-1) ) * dt / ( 2 * dx^2 ) - d(i) * dt / ( 2 * dx^2 ) ...
        - aconst * d(i) * ( xan(i+1) - xan(i-1) ) * dt / ( 2 * R * tempk * dx^2 );
        f_term2(i) = 1 + d(i) * dt / dx^2 ...
        + aconst * ( d(i+1) - d(i-1) ) * ( xan(i+1) - xan(i-1) ) * dt / ( 2 * R * tempk * dx^2) ...
        + aconst * d(i) * ( xan(i+1) - 2 * xan(i) + xan(i-1) ) * dt / ( 2 * R * tempk * dx^2);
        f_term3(i) = - ( d(i+1) - d(i-1) ) * dt / ( 2 * dx^2 ) - d(i) * dt / ( 2 * dx^2 ) ...
        + aconst * d(i) * ( xan(i+1) - xan(i-1) ) * dt / ( 2 * R * tempk * dx^2 );
    end
    # 3重対角行列を作成
    matA = zeros(nx, nx);
    for i = 1:nx
        for j = 1:nx
            if i - j == 1
                matA(i,j) = f_term1(i);
            elseif i == j
                matA(i,j) = f_term2(i);
            elseif i - j == -1
                matA(i,j) = f_term3(i);
            end
        end
    end
    # 2<=i<=nx-1となるように，最初と最後を削除
    matA = cut_mat(matA);
end

# 右辺の行列を計算
function matU = calc_matU(nx, dt, dx, d, xan, tempc, f_term1, f_term3, pmatU, imatU)
    # パラメータを指定
    aconst = -26100; # Bindeman et al. (1998)
    R = 8.314; # gas constant
    tempk = tempc + 273.15;
    # 係数を計算
    c_term1 = zeros(nx);
    c_term2 = zeros(nx);
    c_term3 = zeros(nx);
    for i = 2:nx-1
        c_term1(i) = - ( d(i+1) - d(i-1) ) * dt / ( 2 * dx^2 ) + d(i) * dt / ( 2 * dx^2 ) ...
        + aconst * d(i) * ( xan(i+1) - xan(i-1) ) * dt / ( 2 * R * tempk * dx^2 );
        c_term2(i) = 1 - d(i) * dt / dx^2 ...
        - aconst * ( d(i+1) - d(i-1) ) * ( xan(i+1) - xan(i-1) ) * dt / ( 2 * R * tempk * dx^2) ...
        - aconst * d(i) * ( xan(i+1) - 2 * xan(i) + xan(i-1) ) * dt / ( 2 * R * tempk * dx^2);
        c_term3(i) = ( d(i+1) - d(i-1) ) * dt / ( 2 * dx^2 ) + d(i) * dt / ( 2 * dx^2 ) ...
        - aconst * d(i) * ( xan(i+1) - xan(i-1) ) * dt / ( 2 * R * tempk * dx^2 );
    end
    
    matU = zeros(nx,1);
    for i = 2:nx-1
        if i == 2
            matU(i) = -f_term1(2) * imatU(1) + c_term1(i) * pmatU(i-1) + c_term2(i) * pmatU(i) + c_term3(i) * pmatU(i+1);
        elseif i == nx - 1
            matU(i) = c_term1(i) * pmatU(i-1) + c_term2(i) * pmatU(i) + c_term3(i) * pmatU(i+1) - f_term3(nx-1) * imatU(end);
        else
            matU(i) = c_term1(i) * pmatU(i-1) + c_term2(i) * pmatU(i) + c_term3(i) * pmatU(i+1);
        end
    end
    matU = cut_mat(matU);
end

# メイン処理
function main()
    # CSVファイルの読込
    [distance, xan, imatU] = get_data();

    # パラメータの設定
    tempc = 1000;
    dx = 5*10^(-6);
    dt = 60*60*24;
    nt = 365*300;
    nx = length(distance);
    
    # 拡散係数の計算
    d = diff_coef(xan, tempc, dx, dt);
    [matA, f_term1, f_term3] = calc_matA(nx, dt, dx, d, xan, tempc);

    result = zeros(nx, nt+1);
    result(:,1) = distance;
    
    for j = 1:nt
        if j == 1
            pmatU = imatU;
            matU = calc_matU(nx, dt, dx, d, xan, tempc, f_term1, f_term3, pmatU, imatU);
            fmatU = inv(matA) * matU;
            pmatU = [ imatU(1); fmatU; imatU(end) ];
            result(:,j+1) = pmatU;
        else
            matU = calc_matU(nx, dt, dx, d, xan, tempc, f_term1, f_term3, pmatU, imatU);
            fmatU = inv(matA) * matU;
            pmatU = [ imatU(1); fmatU; imatU(end) ];
            result(:,j+1) = pmatU;
        end
    end
    headline = zeros(1,nt+1);
    result = [ headline; result ];
    csvwrite('tmp.csv', result)

end

# 実行
main()