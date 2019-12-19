#include "commonsettings.h"

namespace Common
{
    //void Settings::set_all()
    //{
        /*fp_tolerance_ = vm_["fp-tolerance"].as<double>();
        hole_cut_ = vm_["hole-cut"].as<std::string>();
        printmesh_ = vm_["printmesh"].as<bool>();
        verbose_ = vm_["verbose"].as<bool>();
        load_estim_ = vm_["load-estim"].as<int>();
        nstripex_ = vm_["nstripex"].as<int>();
        nstripey_ = vm_["nstripey"].as<int>();
        nstripez_ = vm_["nstripez"].as<int>();
        refine_ = vm_["refine"].as<bool>();
        printlm_ = vm_["printlm"].as<bool>();
        print_lm_mesh_ = vm_["print-lm-mesh"].as<bool>();
        print_spc_mesh_ = vm_["print-lm-mesh"].as<bool>();
        refine_tol_ = vm_["refine-tol"].as<double>();
        transfer_cand_donor_ = vm_["transfer-cand-donor"].as<bool>();
        merge_donor_info_ = vm_["merge-donor-info"].as<bool>();
        can_remake_ = vm_["can-remake"].as<bool>();
        remap_batch_ = vm_["remap-batch"].as<bool>();
        force_remake_ = vm_["force-remake"].as<bool>();
        refine_limit_ = vm_["refine-limit"].as<int>();
        greedy_ = vm_["greedy"].as<bool>();
        print_to_file_ = vm_["print-to-file"].as<bool>();
        balance_ = vm_["balance"].as<bool>();
        arearep_ = vm_["area-rep"].as<std::string>();
        dsspat_ = vm_["dsspat"].as<std::string>();
        make_load_balance_ = vm_["make-load-balance"].as<bool>();
        mergebins_ = vm_["mergebins"].as<bool>();*/
    //}

    //void Settings::config(int argc, char* argv[])
    //{
        /*namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description general{"General"};
        general.add_options()
            ("help", "Produce help message")
            ("fp-tolerance", po::value<double>()->default_value(1e-12), "Floating point tolerance")
            ("hole-cut", po::value<std::string>()->default_value("holemap"), "Hole cutting method")
            ("printmesh", po::value<bool>()->default_value(false), "Print meshes")
            ("verbose", po::value<bool>()->default_value(false), "Verbosity")
            ("load-estim", po::value<int>()->default_value(2), "Load estimation method (Default: solver)")
            ("merge-donor-info", po::value<bool>()->default_value(false), "Merge donor info in master")
            ("greedy", po::value<bool>()->default_value(false), "Greedy")
            ("print-to-file", po::value<bool>()->default_value(true), "Generate and write data to files")
            ("balance", po::value<bool>()->default_value(true), "Refine until load balance")
            ("dsspat", po::value<std::string>()->default_value("rm"), "Spatial structure for donor search")
            ("make-load-balance", po::value<bool>()->default_value(true), "Make load balance (for area or hybrid load estimation)")
            ("mergebins", po::value<bool>()->default_value(false), "Merge bins of quad-tree after refinements")
            ;

        po::options_description desc{"Load-map"};
        desc.add_options()
            ("nstripex", po::value<int>()->default_value(2), "Initial number of stripes in x-direction")
            ("nstripey", po::value<int>()->default_value(2), "Initial number of stripes in y-direction")
            ("nstripez", po::value<int>()->default_value(2), "Initial number of stripes in z-direction")
            ("refine", po::value<bool>()->default_value(true), "Can refine load-map")
            ("printlm", po::value<bool>()->default_value(false), "Print load-map")
            ("print-lm-mesh", po::value<bool>()->default_value(false), "Print meshes in load-map")
            ("refine-tol", po::value<double>()->default_value(10.), "Refinement tolerance")
            ("can-remake", po::value<bool>()->default_value(true), "Remake loadmap if degraded")
            ("force-remake", po::value<bool>()->default_value(false), "Force remake of loadmap in each time step (for developer)")
            ("refine-limit", po::value<int>()->default_value(1200), "Maximum number of refinement")
            ("area-rep", po::value<std::string>()->default_value("aabb"), "Area representation in area-based estimation")
            ;

        po::options_description group_spc{"SPC"};
        group_spc.add_options()
            ("print-spc-mesh", po::value<bool>()->default_value(false), "Print meshes in spc")
            ("transfer-cand-donor", po::value<bool>()->default_value(false), "Transfer non-resident candidate donors")
            ("remap-batch", po::value<bool>()->default_value(false), "Make all sp's once when sending in remapping. Batch remapping is faster but tough on memory")
            ;

        all_options.add(general).add(desc).add(group_spc);

        std::ifstream settings_file("settings.ini");

        po::store(po::parse_command_line(argc, argv, all_options), vm_);
        po::store(po::parse_config_file(settings_file, all_options, true), vm_);
        po::notify(vm_);*/
    //}
}
