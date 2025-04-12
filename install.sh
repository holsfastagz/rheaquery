[ -d "$HOME/.local/bin" ] || mkdir $HOME/.local/bin
mkdir $HOME/.local/bin/rheaquery
cp -r src $HOME/.local/bin/rheaquery
cp LICENSE $HOME/.local/bin/rheaquery

echo 'export PATH="$HOME/.local/bin/rheaquery/src:$PATH"' >> $HOME/.bashrc

echo "Rheaquery downloaded successfully. Please run `source ~/.bashrc` for full installation."
