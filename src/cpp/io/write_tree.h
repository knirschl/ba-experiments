//
// Created by knirschl on 25.05.23.
//

#ifndef BA_WRITE_TREE_H
#define BA_WRITE_TREE_H

#include <fstream>
#include "../nj/tree.h"

bool write_newick(Tree const& tree, std::ofstream& writer) {
    if (writer.is_open() && writer.good()) {
        writer << tree.to_newick() << '\n';
    }
    return !writer.fail();
}

bool write_newick(Tree const& tree, std::string const& out_file) {
    std::ofstream streamed_file{out_file};
    bool is_written = write_newick(tree, streamed_file);
    streamed_file.close();
    return is_written;
}

#endif //BA_WRITE_TREE_H
