# del dirs with name
find . -type d -name rf_\* -exec rm -rf {} \;
# del files/dirs with name
find . -name \*fastme_stat\* -delete
# del files matching regex
find . -type f -regex '.*/gene_trees/ba\.exp.[0-9].*\.newick' -delete
# rename files & dirs
find . -name "*.generax_pick*" -exec rename ".generax_pick" "" {} \;
