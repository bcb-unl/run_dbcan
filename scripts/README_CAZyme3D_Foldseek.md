# 使用 CAZyme3D PDB 构建 Foldseek 数据库

你已有：**178373 个 PDB 文件** 和 **cazyme3d_full.tsv**（含 cazyid、family、uniprot_mapped 等列）。按下面步骤生成 run_dbcan 结构搜索所需的 Foldseek 库和映射表。

## 1. 确认 PDB 文件名对应哪一列

- 若 PDB 名为 **RefSeq/CAZyme3D ID**（如 `NP_652145.1.pdb`、`AKV90591.1.pdb`），映射用 **cazyid**。
- 若 PDB 名为 **UniProt**（如 `A0A0G3V1Z8.pdb`），映射用 **uniprot_mapped**。

用你实际目录里的几个文件名看一下即可确定。

## 2. 生成映射表 cazyme3d_to_cazy.tsv

run_dbcan 需要 **target_id → CAZy 家族** 的 TSV（两列：target、CAZy ID），且 **target 必须与 PDB 文件名（去掉 .pdb）一致**。

在仓库根目录或任意位置执行（需安装 pandas）：

```bash
# 若 PDB 以 cazyid 命名（如 NP_652145.1.pdb）
python scripts/cazyme3d_tsv_to_mapping.py /path/to/cazyme3d_full.tsv -o cazyme3d_to_cazy.tsv --id-column cazyid

# 若 PDB 以 UniProt 命名（如 A0A0G3V1Z8.pdb）
python scripts/cazyme3d_tsv_to_mapping.py /path/to/cazyme3d_full.tsv -o cazyme3d_to_cazy.tsv --id-column uniprot_mapped --filter-missing-id
```

`--filter-missing-id` 会去掉 uniprot_mapped 为空（NULL）的行，避免无效 target。

## 3. 用 Foldseek 从 PDB 建库并建索引

在 **同一目录** 下建库（该目录会生成库文件，之后可整体挪到 `db_dir`）：

```bash
# 假设 PDB 所在目录为 /path/to/pdb_dir，共 178373 个 .pdb 文件
PDB_DIR=/path/to/pdb_dir
DB_PREFIX=CAZyme3D

# 3.1 创建 Foldseek 数据库（会生成 CAZyme3D、CAZyme3D_h、CAZyme3D_a3m.ffindex 等）
foldseek createdb "$PDB_DIR" "$DB_PREFIX"

# 3.2 创建索引（加速搜索；大库建议做）
foldseek createindex "$DB_PREFIX" /tmp/foldseek_index --split-memory-limit 40
```

`createindex` 的第二个参数是索引临时目录，可改成你有写权限的路径。若内存紧张，可调小 `--split-memory-limit`（单位 GB）。

## 4. 放入 run_dbcan 的 db_dir

把建库得到的 **所有以 `CAZyme3D` 为前缀的文件** 和 **cazyme3d_to_cazy.tsv** 都放到 run_dbcan 的 `--db_dir` 下，例如：

```text
db_dir/
  CAZyme3D           # createdb 生成
  CAZyme3D_h
  CAZyme3D_a3m.ffindex
  CAZyme3D_a3m.ffdata
  CAZyme3D_dbtype
  ...                # createindex 可能还有其它文件
  cazyme3d_to_cazy.tsv
  CAZy.dmnd          # 其它 run_dbcan 库
  dbCAN.hmm
  ...
```

常量里默认的库前缀是 `CAZyme3D`、映射表名为 `cazyme3d_to_cazy.tsv`，与上面一致即可。

## 5. 运行结构搜索

```bash
run_dbcan CAZyme_annotation \
  --input_raw_data your_proteins.faa \
  --output_dir ./out \
  --db_dir /path/to/db_dir \
  --mode protein \
  --methods diamond,hmm,dbCANsub,structure
```

若只跑结构搜索：`--methods structure`（仍需有 diamond/hmm/dbCANsub 的空结果或同一次跑全方法以生成 overview）。

## 常见问题

- **PDB 名和 TSV 对不上**：务必用 `--id-column cazyid` 或 `--id-column uniprot_mapped` 与真实 PDB 文件名一致；UniProt 时建议加 `--filter-missing-id`。
- **family 列含多域**：如 `CBM50,GH23` 会整段作为 CAZy ID 写入映射表，run_dbcan 会原样使用。
- **库很大**：178k 结构建库和建索引会占不少磁盘与内存，建议在算力节点上执行；createindex 失败时可减小 `--split-memory-limit` 或查 Foldseek 文档。
