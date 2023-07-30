find . -type f -name ba.\* -exec rm {} \;
find . -type f -name rf_\* -exec rm {} \;
find . -name "*.generax_pick*" -exec rename ".generax_pick" "" {} \;
find . -type f -regex '.*/gene_trees/ba\.exp.[0-9].*\.newick' -delete
