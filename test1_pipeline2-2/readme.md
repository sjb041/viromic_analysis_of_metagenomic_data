想法1：分组UViGs聚类成vOTU（(ANI) ≥ 95%+ (AF) ≥ 85%），提取出votu队长， 分组与H的votu代表比较。满足 (ANI ≥ 99.4% over AF ≥ 85%)的Two vOTU队长 were considered clones。【关键点：克隆数基于votu队长的数量】

---

**test1：基于病毒识别的方法`pipeline2-2`**

① 每组样本合并，得到6个文件，代表6组uvigs

```shell
mkdir 01_cat && cd 01_cat

# 获取单独样本的uvigs
cp $HOME/2025_2ME/03_virus_identification/03_pipeline2-2/07_all_sample_results/*_uvigs.fa ./

# 合并每组样本
for i in F1 F2 FG L1 L2 H;do cat ${i}*.fa > ${i}_uvigs.fa;done

# 验证合并前后序列数
grep -c '^>' F1*.fa F2*.fa FG*.fa L1*.fa L2*.fa H*.fa all_sample_uvigs.fa && cd ..
```

---

② 6组uvigs单独聚类（95）

```shell
mkdir -p 02_cluster && cd 02_cluster

for i in F1 F2 FG L1 L2 H
do
    cluster.sh -i ../01_cat/${i}_uvigs.fa -t ../01_cat/${i}_uvigs.fa -a 95 -o ${i}
done

cd ..
```

---

③ 计算克隆（99.4）（克隆数基于==votu代表序列==）

- 方法1（将F的votu代表与H的votu代表合并成一个文件，然后对这个文件99.4聚类）

```shell
mkdir -p 03_clone/method_1 && cd 03_clone/method_1

for i in F1 F2 FG L1 L2
do
	mkdir ${i}-H
	cat ../../02_cluster/${i}/votus.fa ../../02_cluster/H/votus.fa > ${i}-H/cat_votus.fa
    cluster.sh -i ${i}-H/cat_votus.fa -t ${i}-H/cat_votus.fa -a 99.4 -o ${i}-H
done

# 使用聚类结果 cluster.tsv 计算克隆数
for i in F1 F2 FG L1 L2
do 
	Rscript calc_clone-clusters.R -t ../../02_cluster/$i/votus.fa -r ../../02_cluster/H/votus.fa -c ${i}-H/clusters.tsv
done


# 或者使用比对结果 ani.tsv 计算克隆数
for i in F1 F2 FG L1 L2
do
	Rscript calc_clone-ani.R -t ../../02_cluster/$i/votus.fa -r ../../02_cluster/H/votus.fa -a ${i}-H/ani.tsv
done

cd ..
```

---

- 方法2（F的votu代表作为输入序列，H的votu代表作为数据库）

```shell
mkdir -p 03_clone/method_2 && cd 03_clone/method_2

for i in F1 F2 FG L1 L2
do
	mkdir ${i}-H
    cluster.sh -i ../../02_cluster/${i}/votus.fa -t ../../02_cluster/H/votus.fa -a 99.4 -o ${i}-H
done

# 使用比对结果 ani.tsv 计算克隆数
for i in F1 F2 FG L1 L2
do
	Rscript calc_clone-ani.R -t ../../02_cluster/$i/votus.fa -r ../../02_cluster/H/votus.fa -a ${i}-H/ani.tsv
done

cd ..
```

---

- 方法3（H的votu代表作为输入序列，F的votu代表作为数据库）

``` shell
mkdir -p 03_clone/method_3 && cd 03_clone/method_3

for i in F1 F2 FG L1 L2
do
	mkdir ${i}-H
    cluster.sh -i ../../02_cluster/H/votus.fa -t ../../02_cluster/${i}/votus.fa -a 99.4 -o ${i}-H
done

for i in F1 F2 FG L1 L2
do
	Rscript calc_clone-ani.R -t ../../02_cluster/$i/votus.fa -r ../../02_cluster/H/votus.fa -a ${i}-H/ani.tsv
done

cd ..
```

