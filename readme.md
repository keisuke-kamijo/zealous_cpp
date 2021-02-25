# 概要
このリポジトリは[この論文](https://ieeexplore.ieee.org/document/5708146)で提案されたアルゴリズムを実装したものです．
テキストファイルを受け取り，差分プライベートなヒストグラムを出力します．

# コンパイル
c++17でコンパイルできる環境が必要です．

```
make zealous
```

# 実行
実行には4つのコマンドライン引数が必要です．入力ファイル名，εの値，δの値，有用性パラメータを指定してください．有用性パラメータは，削除されるデータを減らしたいときはo,削除されるデータは増えるが頻度の精度を重視したいときはmを入力してください．

例：
```
./zealous input.txt 1 0.001 o
```

#入力ファイル
入力ファイルは(ID,それに属する値)のコンマ区切りで与えてください．右側の値のヒストグラムを作成し出力します．IDはアルゴリズムのステップ1で使用されますが，公開はされません．

