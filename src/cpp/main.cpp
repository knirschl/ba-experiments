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
                   [](auto id) { return std::make_shared<Leaf>(id); });
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
        double div{100.0};
        int step{5};
        for (int i{}; i < 1 * div; i += step) {
            if (i == 1) {
                step = 25;
            }
            matrix_t<double> scaleMatrix{};
            matrix_t<double> sumMatrix{};
            double scale{i / div};
            matscale(species_tree_mat, scale, scaleMatrix);
            matadd(alignment_mat, scaleMatrix, sumMatrix);
            //std::cout << "Scaled Species-Tree-Matrix + Alignment-Matrix =\n" << matstr(sumMatrix) << "\n\n";
            std::shared_ptr<Tree> tree = neighborJoining<>(sumMatrix, alignment_start_leafs);
            std::cout << "Neighbor-joined tree: " << to_newick(*tree) << std::endl;
            // double to string without trailing zeros
            std::ostringstream oss;
            oss << std::setprecision(8) << std::noshowpoint << scale;
            write_newick(*tree, getP(cli_parser) + oss.str() + "S+G.geneTree.newick");
        }
        step = 25;
        for (int i{int(1 * div)}; i <= 2 * div; i += step) {
            matrix_t<double> scaleMatrix{};
            matrix_t<double> sumMatrix{};
            double scale{i / div};
            matscale(species_tree_mat, scale, scaleMatrix);
            matadd(alignment_mat, scaleMatrix, sumMatrix);
            //std::cout << "Scaled Species-Tree-Matrix + Alignment-Matrix =\n" << matstr(distMatrix) << "\n\n";
            std::shared_ptr<Tree> tree = neighborJoining<>(sumMatrix, alignment_start_leafs);
            std::cout << "Neighbor-joined tree: " << to_newick(*tree) << std::endl;
            // double to string without trailing zeros
            std::ostringstream oss;
            oss << std::setprecision(8) << std::noshowpoint << scale;
            write_newick(*tree, getP(cli_parser) + oss.str() + "S+G.geneTree.newick");
        }
        for (int i{int(2.5 * div)}; i <= 10 * div; i *= 2) {
            matrix_t<double> scaleMatrix{};
            matrix_t<double> sumMatrix{};
            double scale{i / div};
            matscale(species_tree_mat, scale, scaleMatrix);
            matadd(alignment_mat, scaleMatrix, sumMatrix);
            //std::cout << "Scaled Species-Tree-Matrix + Alignment-Matrix =\n" << matstr(distMatrix) << "\n\n";
            std::shared_ptr<Tree> tree = neighborJoining<>(sumMatrix, alignment_start_leafs);
            std::cout << "Neighbor-joined tree: " << to_newick(*tree) << std::endl;
            // double to string without trailing zeros
            std::ostringstream oss;
            oss << std::setprecision(8) << std::noshowpoint << scale;
            write_newick(*tree, getP(cli_parser) + oss.str() + "S+G.geneTree.newick");
        }
    }

    return 0;
}


