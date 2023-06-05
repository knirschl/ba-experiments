//
// Created by knirschl on 05.06.23.
//

#include "test.h"

#include <iostream>
#include <memory>
#include "common_types.h"
#include "io/parse_phylip.h"
#include "io/write_tree.h"
#include "nj/matrix.h"
#include "nj/tree.h"
#include "nj/NJSimple.h"

auto convert(vector_t<std::string> string_ids, vector_t<std::shared_ptr<Tree>>& trees) {
    std::transform(string_ids.begin(), string_ids.end(), std::back_inserter(trees),
                   [](auto id) { return std::make_shared<Leaf>(id); });
}

int test() {
    {
        // read species tree
        auto species_tree_pair = parse_from_file<double>(
                "/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_dtl_s20_f100_sites200_GTR_bl1.0_d0.0_l0.0_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed42/species_trees/speciesTree.matrix-sorted.phy");
        auto species_tree_mat = species_tree_pair.first;
        vector_t <std::shared_ptr<Tree>> species_tree_start_leafs;
        convert(species_tree_pair.second, species_tree_start_leafs);
        //std::cout << "Read species tree:\n" << matstr(species_tree_mat) << "\n\n";
        // read alignment
        std::string family_path = "/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_dtl_s20_f100_sites200_GTR_bl1.0_d0.0_l0.0_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed42/families/family_";
        std::string family = "100/";
        std::string alignment_file = "alignment.msa.matrix-sorted.phy";
        std::string alignment_path = family_path + family + alignment_file;
        auto alignment_pair = parse_from_file<double>(alignment_path);
        auto alignment_mat = alignment_pair.first;
        vector_t <std::shared_ptr<Tree>> alignment_start_leafs;
        convert(alignment_pair.second, alignment_start_leafs);
        //std::cout << "Read alignment:\n" << matstr(alignment_mat) << "\n\n";

        // calculate
        matrix_t<double> distMatrix;
        for (double scale{}; scale <= 1.0; scale += 0.1) {
            matscale(species_tree_mat, scale, distMatrix);
            matadd(alignment_mat, distMatrix, distMatrix);
            //std::cout << "Scaled Species-Tree-Matrix + Alignment-Matrix =\n" << matstr(distMatrix) << "\n\n";
            std::shared_ptr <Tree> tree = neighborJoining<>(distMatrix, alignment_start_leafs);
            std::cout << "Neighbor-joined tree: " << to_newick(*tree) << std::endl;
            write_newick(*tree,
                         family_path + family + "gene_trees/ba." + std::to_string(scale) + "S+G.geneTree.newick");
        }
    }
    /*
    {
        // example from de.wiki
        auto dm_id_pair = parse_from_file<double>("/home/fili/Documents/KIT/2023/BA/code/src/cpp/nj/de.wiki.phy");
        auto distMatrix = dm_id_pair.first;
        vector_t<std::shared_ptr<Tree>> trees{};
        std::transform(dm_id_pair.second.begin(), dm_id_pair.second.end(), std::back_inserter(trees),
                       [](auto str_id) { return std::make_shared<Leaf>(str_id); });

        std::shared_ptr<Tree> tree = neighborJoining(distMatrix, trees);
        std::cout << "Neighbor-joined de-tree: " << to_newick(*tree) << std::endl;
    }
     */
    /*
    {
        // example from en.wiki
        auto dm_id_pair = parse_from_file<double>("/home/fili/Documents/KIT/2023/BA/code/src/cpp/nj/en.wiki.phy");
        auto distMatrix = dm_id_pair.first;
        vector_t<std::shared_ptr<Tree>> trees{};
        std::transform(dm_id_pair.second.begin(), dm_id_pair.second.end(), std::back_inserter(trees),
                       [](auto str_id) { return std::make_shared<Leaf>(str_id); });

        std::shared_ptr<Tree> tree = neighborJoining(distMatrix, trees);
        std::cout << "Neighbor-joined en-tree: " << to_newick(*tree) << std::endl;
    }
     */
    return 0;
}