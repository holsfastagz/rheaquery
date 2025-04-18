[ -d "$HOME/.local/bin" ] || mkdir $HOME/.local/bin

[ -d "$HOME/.local/bin/rheaquery" ] && rm -rf "$HOME/.local/bin/rheaquery"
[ -d "$HOME/.local/bin/rheaquery" ] || mkdir $HOME/.local/bin/rheaquery

rm -r * $HOME/.local/bin/rheaquery/*
cp -r src $HOME/.local/bin/rheaquery
cp LICENSE $HOME/.local/bin/rheaquery

grep 'export PATH="$HOME/.local/bin/rheaquery/src:$PATH"' $HOME/.bashrc
|| echo 'export PATH="$HOME/.local/bin/rheaquery/src:$PATH"' >> $HOME/.bashrc

echo "Rheaquery downloaded successfully. Please run \`source ~/.bashrc\` for full installation."
