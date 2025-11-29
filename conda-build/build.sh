#!/bin/bash
set -eu

PKG_DIR="$PREFIX/share/$PKG_NAME-$PKG_VERSION"
BIN_DIR="$PREFIX/bin"
mkdir -p "$PKG_DIR" "$BIN_DIR"

# 1. 拷贝核心文件
for item in bin module main.nf nextflow.config; do
    cp -r "$SRC_DIR/$item" "$PKG_DIR/"
done

# 2. 安装 CRAN / BioC 包（conda-forge 允许）
Rscript -e "install.packages('presto', repos='https://cloud.r-project.org')"
Rscript -e "install.packages('bbknnR', repos='https://cloud.r-project.org')"
Rscript -e "install.packages('openxlsx', repos='https://cloud.r-project.org')"

# 3. 生成入口
cat > "$BIN_DIR/scLine" << 'EOF'
#!/bin/bash
exec nextflow run "$CONDA_PREFIX/share/scline-1.0.0/main.nf" "$@"
EOF
chmod +x "$BIN_DIR/scLine"