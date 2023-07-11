//
// Created by knirschl on 18.04.23.
//

#include <sstream>
#include <iomanip>
#include "io/parse_cli.h"
#include "io/parse_file.h"
#include "io/write_tree.h"
#include "nj/tree.h"
#include "nj/matrix.h"
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

int get_idx(std::vector<std::string> vec, const std::string &key) {
    auto it = find(vec.begin(), vec.end(), key);

    return it != vec.end() ? it - vec.begin() : -1;
}

void run(double scale, const std::shared_ptr<Tree> &tree, std::vector<int> &active,
         const std::vector<std::vector<double>> &species_tree_mat,
         const std::vector<std::basic_string<char>> &species_tree_ids,
         const std::vector<std::vector<double>> &alignment_mat,
         const std::vector<std::basic_string<char>> &alignment_ids,
         argparse::ArgumentParser &cli_parser) {
    auto speciation_pairs{tree->get_speciation_pairs()};
    dist_matrix_t corrected_matrix{alignment_mat};
    for (auto &pair: speciation_pairs) {
        // TODO mapping
        std::string locus1 = idx2leafname[pair.first];
        std::string locus2 = idx2leafname[pair.second];
        int ai = get_idx(alignment_ids, locus1);
        int aj = get_idx(alignment_ids, locus2);
        std::string species1 = leafname2groupname[locus1];
        std::string species2 = leafname2groupname[locus2];
        int si = get_idx(species_tree_ids, species1);
        int sj = get_idx(species_tree_ids, species2);

        corrected_matrix[ai][aj] += scale * species_tree_mat[si][sj];
        corrected_matrix[ai][aj] /= scale + 1;
        corrected_matrix[aj][ai] += scale * species_tree_mat[si][sj];
        corrected_matrix[aj][ai] /= scale + 1;
    }
    neighborJoining<>(corrected_matrix, tree, active);
    std::cout << "Neighbor-joined tree: " << tree->to_newick() << std::endl;

    // double to string without trailing zeros
    std::ostringstream oss;
    oss << std::setprecision(8) << std::noshowpoint << scale;
    write_newick(*tree, getP(cli_parser) + oss.str() + "S~G.geneTree.newick");
}

int main(int argc, char *argv[]) {
    auto cli_parser = build_parser("thesis", "0.1");
    parse(cli_parser, argc, argv);
    //std::cout << getS(cli_parser) << "\n" << getA(cli_parser) << "\n" << getP(cli_parser) << "\n";

    // read species tree
    auto species_tree_pair = parse_phylip_mat_from_file<dist_t>(getS(cli_parser));
    auto species_tree_mat = species_tree_pair.first;
    auto species_tree_ids = species_tree_pair.second;

    // read alignment
    auto alignment_pair = parse_phylip_mat_from_file<dist_t>(getA(cli_parser));
    auto alignment_mat = alignment_pair.first;
    auto alignment_ids = alignment_pair.second;

    auto map_config{getMappingConfig(cli_parser)};

    // calculate
    {
        std::shared_ptr<Tree> tree = reset(species_tree_ids, alignment_ids, map_config);
        std::vector<int> active{leaf_indices};

        // NJ gene tree with only alignment matrix -> old: 0S+G
        neighborJoining<>(alignment_mat, tree, active);
        tree->reroot_APro();
        std::shared_ptr<Tree> backup_tree{tree}; // need to reset tree after each iteration

        // NJ gene tree with corrected values
        double div{100.0};
        int step{5};
        for (int i{}; i < 1 * div; i += step) {
            tree = reset(species_tree_ids, alignment_ids, map_config, backup_tree);
            active = leaf_indices;

            if (i == 1) {
                step = 25;
            }
            run(i / div, tree, active, species_tree_mat, species_tree_ids, alignment_mat,
                alignment_ids, cli_parser);
        }
        step = 25;
        for (int i{int(1 * div)}; i <= 2 * div; i += step) {
            tree = reset(species_tree_ids, alignment_ids, map_config);
            active = leaf_indices;

            run(i / div, tree, active, species_tree_mat, species_tree_ids, alignment_mat,
                alignment_ids, cli_parser);
        }
        for (int i{int(2.5 * div)}; i <= 10 * div; i *= 2) {
            tree = reset(species_tree_ids, alignment_ids, map_config);
            active = leaf_indices;

            run(i / div, tree, active, species_tree_mat, species_tree_ids, alignment_mat,
                alignment_ids, cli_parser);
        }
    }

    return 0;
}


