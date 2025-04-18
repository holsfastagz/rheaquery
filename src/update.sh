echo -e "\nDownloading latest release of Rheaquery. . .\n\n"

# Target directory for download
target_dir="$HOME/.local/bin/rheaquery/downloads"
[ -d "$target_dir" ] || mkdir "$target_dir"

# URL of latest release 
file_url=$(curl -s https://api.github.com/repos/holsfastagz/rheaquery/releases/latest \
           | grep -Eo 'https://[^"]*rheaquery.*\.tar\.gz')

# Filenames/paths
filename=$(basename "$file_url")
filepath="$target_dir/$filename"

# Delete current Rheaquery source
[ -d "$target_dir/rheaquery" ] && rm -rf "$target_dir/rheaquery"

# Download latest Rheaquery release
wget -q -O "$filepath" "$file_url" 

# Unarchive releases and copy them to correct directory
tar -xzf "$filepath" -C "$target_dir"
rm -r "$HOME/.local/bin/rheaquery/src"
cp -r "$target_dir/rheaquery/src" "$HOME/.local/bin/rheaquery"

#rheaquery -v
echo -e "\nRheaquery successfully updated!\n"
