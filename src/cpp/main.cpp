//
// Created by knirschl on 18.04.23.
//

#include <sstream>
#include <iomanip>
#include "io/parse_cli.h"
#include "io/parse_phylip.h"
#include "nj/tree.h"
#include "nj/matrix.h"
#include "nj/NJSimple.h"
#include "io/write_tree.h"
#include <argparse/argparse.hpp>

auto convert(vector_t<std::string> string_ids, vector_t<std::shared_ptr<Tree>>& trees) {
    std::transform(string_ids.begin(), string_ids.end(), std::back_inserter(trees),
                   // TODO mapping.link, this way is very hacky and bad (removes _0_0)
                   [](auto id) { return std::make_shared<Leaf>(id.substr(0, id.length() - 4)); });
}

int main(int argc, char *argv[]) {
    auto cli_parser = build_parser("thesis", "0.1");
    parse(cli_parser, argc, argv);
    //std::cout << getS(cli_parser) << "\n" << getA(cli_parser) << "\n" << getP(cli_parser) << "\n";
    // read species tree
    auto species_tree_pair = parse_from_file<double>(getS(cli_parser));
    auto species_tree_mat = species_tree_pair.first;
    //vector_t<std::shared_ptr<Tree>> species_tree_start_leafs;
    //convert(species_tree_pair.second, species_tree_start_leafs);
    // read alignment
    auto alignment_pair = parse_from_file<double>(getA(cli_parser));
    auto alignment_mat = alignment_pair.first;
    vector_t<std::shared_ptr<Tree>> alignment_start_leafs;
    convert(alignment_pair.second, alignment_start_leafs);

    // calculate

    {
        matrix_t<double> scaleMatrix{};
        matrix_t<double> sumMatrix{};
        // 0S+1G
        std::shared_ptr<Tree> tree = neighborJoining<>(alignment_mat, alignment_start_leafs);
        write_newick(*tree, getP(cli_parser) + "0S+" + "1G.geneTree.newick");
    }
    {
        matrix_t<double> scaleMatrix{};
        matrix_t<double> sumMatrix{};
        // 0.3S+1G
        matscale(species_tree_mat, 0.3, scaleMatrix);
        matadd(scaleMatrix, alignment_mat, sumMatrix);
        std::shared_ptr<Tree> tree = neighborJoining<>(sumMatrix, alignment_start_leafs);
        write_newick(*tree, getP(cli_parser) + "0.3S+" + "1G.geneTree.newick");
    }
    {
        matrix_t<double> scaleMatrix{};
        matrix_t<double> sumMatrix{};
        // 1S+1G
        matadd(species_tree_mat, alignment_mat, sumMatrix);
        std::shared_ptr<Tree> tree = neighborJoining<>(sumMatrix, alignment_start_leafs);
        write_newick(*tree, getP(cli_parser) + "1S+" + "1G.geneTree.newick");
    }
    {
        matrix_t<double> scaleMatrix{};
        matrix_t<double> sumMatrix{};
        // 10S+1G
        matscale(species_tree_mat, 10, scaleMatrix);
        matadd(scaleMatrix, alignment_mat, sumMatrix);
        std::shared_ptr<Tree> tree = neighborJoining<>(sumMatrix, alignment_start_leafs);
        write_newick(*tree, getP(cli_parser) + "10S+" + "1G.geneTree.newick");
    }
    {
        matrix_t<double> scaleMatrix{};
        matrix_t<double> sumMatrix{};
        // 100S+1G
        matscale(species_tree_mat, 100, scaleMatrix);
        matadd(scaleMatrix, alignment_mat, sumMatrix);
        std::shared_ptr<Tree> tree = neighborJoining<>(sumMatrix, alignment_start_leafs);
        write_newick(*tree, getP(cli_parser) + "100S+" + "1G.geneTree.newick");
    }

    /*
    for (double scale{}; scale <= 1.0; scale += 0.1) {
        matscale(species_tree_mat, scale, distMatrix);
        matadd(alignment_mat, distMatrix, distMatrix);
        //std::cout << "Scaled Species-Tree-Matrix + Alignment-Matrix =\n" << matstr(distMatrix) << "\n\n";
        std::shared_ptr <Tree> tree = neighborJoining<>(distMatrix, alignment_start_leafs);
        std::cout << "Neighbor-joined tree: " << to_newick(*tree) << std::endl;
        // double to string without trailing zeros
        std::ostringstream oss;
        oss << std::setprecision(8) << std::noshowpoint << scale;
        write_newick(*tree, getP(cli_parser) + oss.str() + "S+G.geneTree.newick");
    }

    for (int scale{5}; scale <= 100; scale += 5) {
        matscale(species_tree_mat, scale, distMatrix);
        matadd(alignment_mat, distMatrix, distMatrix);
        //std::cout << "Scaled Species-Tree-Matrix + Alignment-Matrix =\n" << matstr(distMatrix) << "\n\n";
        std::shared_ptr <Tree> tree = neighborJoining<>(distMatrix, alignment_start_leafs);
        std::cout << "Neighbor-joined tree: " << to_newick(*tree) << std::endl;
        write_newick(*tree, getP(cli_parser) + std::to_string(scale) + "S+G.geneTree.newick");
    }
    */

    return 0;
}


