//
// Created by knirschl on 18.04.23.
//

#include <iostream>
#include <memory>
#include "common_types.h"
#include "io/parse_phylip.h"
#include "nj/matrix.h"
#include "nj/tree.h"
#include "nj/NJSimple.h"

auto convert(vector_t<std::string> string_ids, vector_t<std::shared_ptr<NTree>>& trees) { 
    std::transform(string_ids.begin(), string_ids.end(), std::back_inserter(trees),
                   [](auto id) { return std::make_shared<NLeaf>(id); });
}

int main() {
   {
        // read species tree
        auto species_tree_pair = parse_from_file<double>("/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_dtl_s20_f100_sites200_GTR_bl1.0_d0.0_l0.0_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed42/species_trees/speciesTree.matrix.phy");
        auto species_tree_mat = species_tree_pair.first;
        vector_t<std::shared_ptr<NTree>> species_tree_start_leafs;
        convert(species_tree_pair.second, species_tree_start_leafs);
        std::cout << "Read species tree:\n" << matstr(species_tree_mat) << "\n\n";
        // read alignment
        std::string family_path = "/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_dtl_s20_f100_sites200_GTR_bl1.0_d0.0_l0.0_t1.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed42/families/family_";
        std::string alignment_file = "alignment.msa.matrix.phy";
        std::string alignment_path = family_path + "100/" + alignment_file;
        auto alignment_pair = parse_from_file<double>(alignment_path);
        auto alignment_mat = alignment_pair.first;
        vector_t<std::shared_ptr<NTree>> alignment_start_leafs;
        convert(alignment_pair.second, alignment_start_leafs);
        std::cout << "Read alignment:\n" << matstr(alignment_mat) << "\n\n";
        // calculate
        matrix_t<double> distMatrix;
        for (double scale{}; scale <= 1.0; scale += 0.1) {
            matscale(species_tree_mat, scale, distMatrix);
            matadd(alignment_mat, distMatrix, distMatrix);
            std::cout << "Calculated distance matrix:\n" << matstr(distMatrix) << "\n\n";
            std::shared_ptr<NTree> tree = neighborJoining<>(distMatrix, alignment_start_leafs);
            std::cout << "Neighbor-joined tree: " << to_fasta(*tree) << std::endl;
        }
    }
    /*
    {
        // example from de.wiki
        auto dm_id_pair = parse_from_file<double>("/home/fili/Documents/KIT/2023/BA/code/src/cpp/nj/de.wiki.phy");
        auto distMatrix = dm_id_pair.first;
        vector_t<std::shared_ptr<NTree>> trees{};
        std::transform(dm_id_pair.second.begin(), dm_id_pair.second.end(), std::back_inserter(trees),
                       [](auto str_id) { return std::make_shared<NLeaf>(str_id); });

        std::shared_ptr<NTree> tree = neighborJoining(distMatrix, trees);
        std::cout << "Neighbor-joined de-tree: " << to_fasta(*tree) << std::endl;
    }
     */
    /*
    {
        // example from en.wiki
        auto dm_id_pair = parse_from_file<double>("/home/fili/Documents/KIT/2023/BA/code/src/cpp/nj/en.wiki.phy");
        auto distMatrix = dm_id_pair.first;
        vector_t<std::shared_ptr<NTree>> trees{};
        std::transform(dm_id_pair.second.begin(), dm_id_pair.second.end(), std::back_inserter(trees),
                       [](auto str_id) { return std::make_shared<NLeaf>(str_id); });

        std::shared_ptr<NTree> tree = neighborJoining(distMatrix, trees);
        std::cout << "Neighbor-joined en-tree: " << to_fasta(*tree) << std::endl;
    }
     */
    return 0;
}
