# scLine/conda-build/build.sh
#!/bin/bash

# 创建目标目录
mkdir -p $PREFIX/share/$PKG_NAME-$PKG_VERSION
mkdir -p $PREFIX/bin

# 拷贝项目文件
cp -r $SRC_DIR/* $PREFIX/share/$PKG_NAME-$PKG_VERSION/

# 创建可执行脚本
cat > $PREFIX/bin/scLine << 'EOF'
#!/bin/bash
exec nextflow run "$CONDA_PREFIX/share/scline-1.0.0/main.nf" "$@"
EOF

chmod +x $PREFIX/bin/scLine