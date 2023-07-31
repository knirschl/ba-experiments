//
// Created by knirschl on 18.04.23.
//

#include <sstream>
#include <iomanip>
#include "io/parse_cli.h"
#include "io/parse_file.h"
#include "io/write_file.h"
#include "nj/tree.h"
#include "nj/NJSimple.h"
#include "misc/meta.h"

/**
 * Resets tree and mappings that need to be reset.
 *
 * ==========
 * INFO: Only needed for figuring out what settings are best
 * ==========
 *
 * @param map_config
 * @return
 */
auto reset(const std::vector<std::string> &species_tree_ids,
           const std::vector<std::string> &alignment_ids,
           const std::tuple<std::string, bool, std::string> &map_config,
           const std::shared_ptr<Tree> &old_tree = nullptr) {
    // fill map
    if (get<0>(map_config).empty()) {
        // no mapping provided
        for (int i{}; i < alignment_ids.size(); i++) {
            leafname2groupname.emplace(alignment_ids[i], species_tree_ids[i]);
        }
    } else {
        leafname2groupname = parse_mapping_from_cfg(map_config);
    }
    // create tree
    std::shared_ptr<Tree> tree = std::make_shared<Tree>();
    tree->make_leafs(alignment_ids);
    tree = old_tree != nullptr ? old_tree : tree;
    return tree;
}

std::string to_string(const vector_t<double> &vec) {
    return "[" + accumulate(vec.begin() + 1, vec.end(),
                            std::to_string(vec[0]), [](const std::string &acc, double b) {
        return acc + ", " + std::to_string(b);
    }) + "]";
}

std::string to_string(const matrix_t<double> &mat) {
    return "[" + accumulate(mat.begin() + 1, mat.end(),
                            to_string(mat[0]), [](const std::string &acc, vector_t<double> b) {
        return acc + ",\n" + to_string(b);
    }) + "]";
}

bool run(double scale, const std::shared_ptr<Tree> &tree, std::vector<int> &active,
         const dist_matrix_t &species_tree_mat,
         const dist_matrix_t &alignment_mat,
         const vector_t<std::string> &alignment_ids,
         argparse::ArgumentParser &cli_parser) {

    //std::cout << "Start tree: " << tree->to_string()
    //          << "\nSpecies tree matrix:\n" << to_string(species_tree_mat)
    //        << "\nAlignment matrix:\n" << to_string(alignment_mat) << "\n";

    dist_matrix_t corrected_matrix{alignment_mat};
    for (auto &pair: tree->get_speciation_pairs()) {
        std::string locus1 = idx2leafname[pair.first];
        std::string locus2 = idx2leafname[pair.second];
        int ai = leafname2matidx[locus1];
        int aj = leafname2matidx[locus2];
        int si = groupname2matidx[leafname2groupname[locus1]];
        int sj = groupname2matidx[leafname2groupname[locus2]];

        corrected_matrix[ai][aj] += scale * species_tree_mat[si][sj];
        corrected_matrix[ai][aj] /= scale + 1;
        corrected_matrix[aj][ai] = corrected_matrix[ai][aj];
    }
    // double to string without trailing zeros
    std::ostringstream oss;
    oss << std::setprecision(8) << std::noshowpoint << scale;
    // matrix or tree
    bool ret{true};
    if (int c{get_c(cli_parser)}) {
        ret = write_phylip(corrected_matrix, alignment_ids, get_output_prefix(cli_parser) + oss.str() + "S~G.matrix.phy");
        if (c == 2) {
            return ret;
        }
    }
    //std::cout << "\nCorrected matrix:\n" << to_string(corrected_matrix) << "\n";
    neighborJoining<>(corrected_matrix, tree, active);

    std::cout << "Neighbor-joined tree (" << oss.str() << "S~G): " << tree->to_newick()
              << std::endl;
    return ret & write_newick(*tree, get_output_prefix(cli_parser) + oss.str() + "S~G.geneTree.newick");
}

/*
-s
/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_DL_s100_f100_sites200_GTR_bl1.0_d1.0_l1.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed1650734/species_trees/speciesTree.matrix.phy
-a
/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_DL_s100_f100_sites200_GTR_bl1.0_d1.0_l1.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed1650734/families/family_042/alignment.msa.matrix.phy
-p
/home/fili/Documents/KIT/2023/BA/code/src/test/output/ba.exp.
-m
/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_DL_s100_f100_sites200_GTR_bl1.0_d1.0_l1.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed1650734/families/family_042/mappings/mapping.link
 */
int main(int argc, char *argv[]) {
    // --- parse ---
    // get cli inputs
    auto cli_parser = build_parser("thesis", "0.2");
    parse(cli_parser, argc, argv);
    //std::cout << getS(cli_parser) << "\n" << getA(cli_parser) << "\n" << getP(cli_parser) << "\n";

    // read species tree
    auto species_tree_pair = parse_phylip_mat_from_file<dist_t>(get_species_matrix(cli_parser));
    auto species_tree_mat = species_tree_pair.first;
    auto species_tree_ids = species_tree_pair.second;
    int s_cnt{};
    for_each(species_tree_ids.begin(), species_tree_ids.end(),
             [&s_cnt](auto &s) { groupname2matidx.emplace(s, s_cnt++); });

    // read alignment
    auto alignment_pair = parse_phylip_mat_from_file<dist_t>(get_alignment_matrix(cli_parser));
    auto alignment_mat = alignment_pair.first;
    auto alignment_ids = alignment_pair.second;
    int a_cnt{};
    for_each(alignment_ids.begin(), alignment_ids.end(),
             [&a_cnt](auto &s) { leafname2matidx.emplace(s, a_cnt++); });

    // read mapping
    auto map_config{get_mapping_config(cli_parser)};

    // calculate
    {
        std::shared_ptr<Tree> tree = reset(species_tree_ids, alignment_ids, map_config);
        active = leaf_indices;
        // NJ gene tree with only alignment matrix (0S+G)
        neighborJoining<>(alignment_mat, tree, active);
    }
    //std::cout << "Start tree : " << tree->to_newick() << "\n" << tree->node_info() << "\n";
    tree->reroot_APro();
    std::shared_ptr<Tree> backup_tree{tree}; // need to reset tree after each iteration
    //std::cout << "A-Pro tree: " << tree->to_newick() << "\n" << tree->node_info() << "\n";

    // --- calculate ---
    // NJ gene tree with corrected values
    double scales[] = { /*
            0, 0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,1.9, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4,
            2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.85, 2.95, 3.1, 3.2, 3.3, 3.4, 3.5, 4, 4.5, 5, 7.5, 10
            */
            0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85,
            1.95, 2.8, 2.9, 3.0, 3.05, 3.15, 3.25, 3.35, 3.45, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.05,
            4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95
            };
    for (double scale: scales) {
        tree = reset(species_tree_ids, alignment_ids, map_config);
        active = leaf_indices;

        run(scale, tree, active, species_tree_mat, alignment_mat, alignment_ids, cli_parser);
    }

    return 0;
}


