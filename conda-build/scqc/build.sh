# #!/usr/bin/env bash
# set -euo pipefail

# # 1. 把 nextflow 主流程放到 $PREFIX/share/scqc/
# mkdir -p $PREFIX/share/scqc
# cp $SRC_DIR/scqc_main.nf   $PREFIX/share/scqc/

# # 2. 生成可执行入口脚本 $PREFIX/bin/scqc
# mkdir -p $PREFIX/bin
# cat > $PREFIX/bin/scqc <<'EOF'
# #!/usr/bin/env bash
# # 自动找到 conda 环境下的脚本
# NXF_HOME=${CONDA_PREFIX}/share/scqc
# exec nextflow "$NXF_HOME/scqc_main.nf" "$@"
# EOF
# chmod +x $PREFIX/bin/scqc