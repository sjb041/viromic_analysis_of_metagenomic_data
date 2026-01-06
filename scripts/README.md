```bash
# 给目录 scripts 下的所有文件加上可执行权限
find $HOME/2025_2ME/scripts -type f | xargs chmod +x

# 把 scripts 和子目录 VMGCPipelines 加入到 PATH 中
echo 'export PATH="$HOME/2025_2ME/scripts:$HOME/2025_2ME/scripts/VMGCPipelines:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
