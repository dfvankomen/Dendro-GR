//
// Created by milinda on 8/23/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief
 */
//

#include "parameters.h"

#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

namespace bssn {

mem::memory_pool<double> BSSN_MEM_POOL = mem::memory_pool<double>(0, 16);
unsigned int BSSN_ELE_ORDER            = 6;
unsigned int BSSN_PADDING_WIDTH        = BSSN_ELE_ORDER >> 1u;

unsigned int BSSN_IO_OUTPUT_FREQ       = 10;
unsigned int BSSN_TIME_STEP_OUTPUT_FREQ =
    std::numeric_limits<unsigned int>::max();

double BSSN_BH_MERGE_TIME          = std::numeric_limits<double>::max();
unsigned int BSSN_BH_MERGE_STEP    = std::numeric_limits<unsigned int>::max();

unsigned int BSSN_GW_EXTRACT_FREQ  = std::numeric_limits<unsigned int>::max();

unsigned int BSSN_REMESH_TEST_FREQ = 10;

unsigned int BSSN_REMESH_TEST_FREQ_AFTER_MERGER        = 10;
unsigned int BSSN_GW_EXTRACT_FREQ_AFTER_MERGER         = 10;
double BSSN_IO_OUTPUT_GAP                              = 1.0;

unsigned int BSSN_USE_WAVELET_TOL_FUNCTION             = 0;
double BSSN_WAVELET_TOL                                = 0.0001;
double BSSN_GW_REFINE_WTOL                             = 1e-4;
double BSSN_WAVELET_TOL_MAX                            = 0.001;
double BSSN_WAVELET_TOL_FUNCTION_R0                    = 10.0;
double BSSN_WAVELET_TOL_FUNCTION_R1                    = 50.0;

double BSSN_CFL_FACTOR                                 = 0.1;

bool BSSN_KO_SIGMA_SCALE_BY_CONFORMAL                  = false;
bool BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY = false;
double BSSN_EPSILON_CAKO_GAUGE                         = 0.99;
double BSSN_EPSILON_CAKO_OTHER                         = 0.3;
bool BSSN_CAKO_ENABLED                                 = false;
double BSSN_CAHD_C                                     = 0.0;

double BSSN_LOAD_IMB_TOL                               = 0.1;
unsigned int BSSN_SPLIT_FIX                            = 2;
unsigned int BSSN_ASYNC_COMM_K                         = 4;
double BSSN_RK_TIME_BEGIN                              = 0;
double BSSN_RK_TIME_END                                = 10;

double BSSN_RK45_DESIRED_TOL                           = 1e-6;

unsigned int BSSN_RK_TYPE;

unsigned int BSSN_CHECKPT_FREQ            = 10;
unsigned int BSSN_RESTORE_SOLVER          = 0;
unsigned int BSSN_ENABLE_BLOCK_ADAPTIVITY = 0;

double BSSN_ETA_R0                        = 1.31;
double BSSN_ETA_POWER[]                   = {2.0, 2.0};

std::string BSSN_VTU_FILE_PREFIX          = "bssn_gr";
std::string BSSN_CHKPT_FILE_PREFIX        = "bssn_cp";
std::string BSSN_PROFILE_FILE_PREFIX      = "bssn_prof";

double BSSN_BH1_AMR_R                     = 2.0;
double BSSN_BH2_AMR_R                     = 2.0;
// ratio for the near to far portion, was originally 2.5
double BSSN_AMR_R_RATIO                   = 2.0;

double BSSN_BH1_CONSTRAINT_R              = 5.0;
double BSSN_BH2_CONSTRAINT_R              = 5.0;

double BSSN_BH1_MASS;
double BSSN_BH2_MASS;

unsigned int BSSN_BH1_MAX_LEV    = std::numeric_limits<unsigned int>::max();
unsigned int BSSN_BH2_MAX_LEV    = std::numeric_limits<unsigned int>::max();

unsigned int BSSN_INIT_GRID_ITER = 10;

BH BH1;
BH BH2;
Point BSSN_BH_LOC[2];
unsigned int BSSN_DIM      = 3;
unsigned int BSSN_MAXDEPTH = 8;
unsigned int BSSN_MINDEPTH = 3;

unsigned int BSSN_ID_TYPE  = 0;

double BSSN_GRID_MIN_X     = -50.0;
double BSSN_GRID_MAX_X     = 50.0;
double BSSN_GRID_MIN_Y     = -50.0;
double BSSN_GRID_MAX_Y     = 50.0;
double BSSN_GRID_MIN_Z     = -50.0;
double BSSN_GRID_MAX_Z     = 50.0;

double BSSN_BLK_MIN_X      = -6.0;
double BSSN_BLK_MIN_Y      = -6.0;
double BSSN_BLK_MIN_Z      = -6.0;

double BSSN_BLK_MAX_X      = 6.0;
double BSSN_BLK_MAX_Y      = 6.0;
double BSSN_BLK_MAX_Z      = 6.0;

double BSSN_COMPD_MIN[3]  = {BSSN_GRID_MIN_X, BSSN_GRID_MIN_Y, BSSN_GRID_MIN_Z};
double BSSN_COMPD_MAX[3]  = {BSSN_GRID_MAX_X, BSSN_GRID_MAX_Y, BSSN_GRID_MAX_Z};

double BSSN_OCTREE_MIN[3] = {0.0, 0.0, 0.0};
double BSSN_OCTREE_MAX[3] = {(double)(1u << BSSN_MAXDEPTH),
                             (double)(1u << BSSN_MAXDEPTH),
                             (double)(1u << BSSN_MAXDEPTH)};

//@note assumes the computational domain is a cube as well.
double BSSN_RK45_TIME_STEP_SIZE = BSSN_CFL_FACTOR *
                                  (BSSN_COMPD_MAX[0] - BSSN_COMPD_MIN[0]) *
                                  (1.0 / (double)(1u << BSSN_MAXDEPTH));

// calculate the minimum dx
double BSSN_CURRENT_MIN_DX = (BSSN_COMPD_MAX[0] - BSSN_COMPD_MIN[0]) *
                             (1.0 / (double)(1u << BSSN_MAXDEPTH));

unsigned int BSSN_LAMBDA[4]                              = {1, 1, 1, 1};
double BSSN_LAMBDA_F[2]                                  = {1.0, 0.0};
double BSSN_TRK0                                         = 0.0;
double ETA_CONST                                         = 2.0;
double ETA_R0                                            = 50.0;
double ETA_DAMPING                                       = 1.0;
double ETA_DAMPING_EXP                                   = 1.0;
double CHI_FLOOR                                         = 0.1;
double KO_DISS_SIGMA                                     = 0.01;

unsigned int RIT_ETA_FUNCTION                            = 1;
double RIT_ETA_OUTER                                     = 0.25;
double RIT_ETA_CENTRAL                                   = 2.0;
double RIT_ETA_WIDTH                                     = 40.0;

unsigned int DISSIPATION_TYPE                            = 0;

unsigned int BSSN_DENDRO_GRAIN_SZ                        = 1000;

double BSSN_DENDRO_AMR_FAC                               = 0.1;
double BSSN_DENDRO_AMR_FAC_POST_MERGER                   = 0.0;

unsigned int BSSN_NUM_REFINE_VARS                        = 1;
unsigned int BSSN_REFINE_VARIABLE_INDICES[BSSN_NUM_VARS] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

unsigned int BSSN_NUM_EVOL_VARS_VTU_OUTPUT               = 1;
unsigned int BSSN_NUM_CONST_VARS_VTU_OUTPUT              = 1;
unsigned int BSSN_VTU_OUTPUT_EVOL_INDICES[BSSN_NUM_VARS] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
unsigned int BSSN_VTU_OUTPUT_CONST_INDICES[BSSN_CONSTRAINT_NUM_VARS] = {
    0, 1, 2, 3, 4, 5};

unsigned int BSSN_XI[3]                         = {0, 0, 0};

unsigned int BSSN_DISSIPATION_NC                = 0;

unsigned int BSSN_DISSIPATION_S                 = 10;

double BSSN_EH_REFINE_VAL                       = 0.3;
double BSSN_EH_COARSEN_VAL                      = 0.4;

// by default use WAMR refinement.
RefinementMode BSSN_REFINEMENT_MODE             = RefinementMode::WAMR;

bool BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE = false;

bool BSSN_VTU_Z_SLICE_ONLY                      = true;

unsigned int BSSN_LTS_TS_OFFSET                 = 4;

bool BSSN_MERGED_CHKPT_WRITTEN                  = false;

double BSSN_CURRENT_RK_COORD_TIME               = 0;
unsigned int BSSN_CURRENT_RK_STEP               = 0;

unsigned int BSSN_NYQUIST_M                     = 0;

bool BSSN_SCALE_VTU_AND_GW_EXTRACTION           = false;

unsigned int BSSN_GW_EXTRACT_FREQ_TRUE          = 0;

unsigned int BSSN_IO_OUTPUT_FREQ_TRUE           = 0;

double BSSN_SSL_SIGMA                           = 20.0;
double BSSN_SSL_H                               = 0.6;

/***@brief: derivs workspace*/
double* BSSN_DERIV_WORKSPACE                    = nullptr;

struct ParameterInformation {
    std::string key;
    std::any value_reference;
    bool required;

    template <typename T>
    ParameterInformation(std::string name, T& var, bool required)
        : key(std::move(name)),
          value_reference(std::ref(var)),
          required(required) {}

    // vector type
    template <typename T>
    ParameterInformation(std::string name, std::vector<T>& var, bool required)
        : key(std::move(name)),
          value_reference(std::ref(var)),
          required(required) {}
};

void readParamTOMLFile(const char* fName, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    auto parFile = toml::parse(fName);
    std::unordered_set<std::string> used_params;

    auto set_param = [&](auto& pardata, const ParameterInformation& param) {
        if (param.required && !pardata.contains(param.key)) {
            throw std::runtime_error("Missing required parameter: " +
                                     param.key);
        }
        // so, if it has the key then we set it, if not we print a warning
        if (pardata.contains(param.key)) {
            // insert the key as used:
            used_params.insert(param.key);
            if (auto* p =
                    std::any_cast<std::reference_wrapper<std::vector<int>>>(
                        &param.value_reference)) {
                p->get() = toml::find<std::vector<int>>(pardata, param.key);
            } else if (auto* p = std::any_cast<
                           std::reference_wrapper<std::vector<unsigned int>>>(
                           &param.value_reference)) {
                std::vector<int> vals_tmp;
                vals_tmp = toml::find<std::vector<int>>(pardata, param.key);
                std::vector<unsigned int> vals_tmp_convert(vals_tmp.begin(),
                                                           vals_tmp.end());

                p->get() = vals_tmp_convert;
            } else if (auto* p = std::any_cast<
                           std::reference_wrapper<std::vector<double>>>(
                           &param.value_reference)) {
                p->get() = toml::find<std::vector<double>>(pardata, param.key);
            } else if (auto* p = std::any_cast<std::reference_wrapper<int>>(
                           &param.value_reference)) {
                p->get() = pardata[param.key].as_integer();
            } else if (auto* p =
                           std::any_cast<std::reference_wrapper<unsigned int>>(
                               &param.value_reference)) {
                p->get() = pardata[param.key].as_integer();
            } else if (auto* p =
                           std::any_cast<std::reference_wrapper<std::string>>(
                               &param.value_reference)) {
                p->get() = pardata[param.key].as_string();
            } else if (auto* p = std::any_cast<std::reference_wrapper<double>>(
                           &param.value_reference)) {
                p->get() = pardata[param.key].as_floating();
            } else if (auto* p = std::any_cast<std::reference_wrapper<float>>(
                           &param.value_reference)) {
                p->get() = pardata[param.key].as_floating();
            } else if (auto* p = std::any_cast<std::reference_wrapper<bool>>(
                           &param.value_reference)) {
                p->get() = pardata[param.key].as_boolean();
            } else {
                throw std::runtime_error("Unsupported parameter type: " +
                                         param.key);
            }
        } else {
            if (rank == 0) {
                std::cout << YLW << "\tOPTIONAL PARAMETER '" << param.key
                          << "' wasn't found! Using default!" << NRM
                          << std::endl;
            }
        }
    };

    std::vector<ParameterInformation> parsList = {
        {"BSSN_IO_OUTPUT_FREQ", bssn::BSSN_IO_OUTPUT_FREQ, true},
        {"BSSN_REMESH_TEST_FREQ", bssn::BSSN_REMESH_TEST_FREQ, true},
        {"BSSN_CHECKPT_FREQ", bssn::BSSN_CHECKPT_FREQ, true},
        // {"BSSN_IO_OUTPUT_GAP", bssn::BSSN_IO_OUTPUT_GAP, true}, // this is
        // actually unused!
        {"BSSN_VTU_FILE_PREFIX", bssn::BSSN_VTU_FILE_PREFIX, true},
        {"BSSN_CHKPT_FILE_PREFIX", bssn::BSSN_CHKPT_FILE_PREFIX, true},
        {"BSSN_PROFILE_FILE_PREFIX", bssn::BSSN_PROFILE_FILE_PREFIX, true},
        {"BSSN_RESTORE_SOLVER", bssn::BSSN_RESTORE_SOLVER, true},
        {"BSSN_ID_TYPE", bssn::BSSN_ID_TYPE, true},
        {"BSSN_ENABLE_BLOCK_ADAPTIVITY", bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY,
         true},
        {"BSSN_BLK_MIN_X", bssn::BSSN_BLK_MIN_X, true},
        {"BSSN_BLK_MIN_Y", bssn::BSSN_BLK_MIN_Y, true},
        {"BSSN_BLK_MIN_Z", bssn::BSSN_BLK_MIN_Z, true},
        {"BSSN_BLK_MAX_X", bssn::BSSN_BLK_MAX_X, true},
        {"BSSN_BLK_MAX_Y", bssn::BSSN_BLK_MAX_Y, true},
        {"BSSN_BLK_MAX_Z", bssn::BSSN_BLK_MAX_Z, true},
        {"BSSN_DENDRO_GRAIN_SZ", bssn::BSSN_DENDRO_GRAIN_SZ, true},
        {"BSSN_ASYNC_COMM_K", bssn::BSSN_ASYNC_COMM_K, true},
        {"BSSN_DENDRO_AMR_FAC", bssn::BSSN_DENDRO_AMR_FAC, true},
        {"BSSN_LOAD_IMB_TOL", bssn::BSSN_LOAD_IMB_TOL, true},
        {"BSSN_RK_TIME_BEGIN", bssn::BSSN_RK_TIME_BEGIN, true},
        {"BSSN_RK_TIME_END", bssn::BSSN_RK_TIME_END, true},
        {"BSSN_RK_TYPE", bssn::BSSN_RK_TYPE, true},
        {"BSSN_RK45_TIME_STEP_SIZE", bssn::BSSN_RK45_TIME_STEP_SIZE, true},
        {"BSSN_RK45_DESIRED_TOL", bssn::BSSN_RK45_DESIRED_TOL, true},
        {"BSSN_DIM", bssn::BSSN_DIM, true},
        {"BSSN_MAXDEPTH", bssn::BSSN_MAXDEPTH, true},
        {"BSSN_GRID_MIN_X", bssn::BSSN_GRID_MIN_X, true},
        {"BSSN_GRID_MIN_Y", bssn::BSSN_GRID_MIN_Y, true},
        {"BSSN_GRID_MIN_Z", bssn::BSSN_GRID_MIN_Z, true},
        {"BSSN_GRID_MAX_X", bssn::BSSN_GRID_MAX_X, true},
        {"BSSN_GRID_MAX_Y", bssn::BSSN_GRID_MAX_Y, true},
        {"BSSN_GRID_MAX_Z", bssn::BSSN_GRID_MAX_Z, true},
        {"ETA_CONST", bssn::ETA_CONST, true},
        {"ETA_R0", bssn::ETA_R0, true},
        {"ETA_DAMPING", bssn::ETA_DAMPING, true},
        {"ETA_DAMPING_EXP", bssn::ETA_DAMPING_EXP, true},
        {"RIT_ETA_FUNCTION", bssn::RIT_ETA_FUNCTION, false},
        {"RIT_ETA_OUTER", bssn::RIT_ETA_OUTER, false},
        {"RIT_ETA_CENTRAL", bssn::RIT_ETA_CENTRAL, false},
        {"RIT_ETA_WIDTH", bssn::RIT_ETA_WIDTH, false},
        {"BSSN_KO_SIGMA_SCALE_BY_CONFORMAL",
         bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL, false},
        {"BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY",
         bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY, false},
        {"BSSN_EPSILON_CAKO_GAUGE", bssn::BSSN_EPSILON_CAKO_GAUGE, false},
        {"BSSN_EPSILON_CAKO_OTHER", bssn::BSSN_EPSILON_CAKO_OTHER, false},
        {"BSSN_CAHD_C", bssn::BSSN_CAHD_C, false},
        // lambda and xi are added later
        {"BSSN_ELE_ORDER", bssn::BSSN_ELE_ORDER, false},
        {"CHI_FLOOR", bssn::CHI_FLOOR, true},
        {"BSSN_TRK0", bssn::BSSN_TRK0, true},
        {"DISSIPATION_TYPE", bssn::DISSIPATION_TYPE, false},
        {"KO_DISS_SIGMA", bssn::KO_DISS_SIGMA, true},
        {"BSSN_ETA_R0", bssn::BSSN_ETA_R0, true},
        // eta power 0 is added later
        {"BSSN_USE_WAVELET_TOL_FUNCTION", bssn::BSSN_USE_WAVELET_TOL_FUNCTION,
         true},
        {"BSSN_WAVELET_TOL", bssn::BSSN_WAVELET_TOL, true},
        {"BSSN_WAVELET_TOL_MAX", bssn::BSSN_WAVELET_TOL_MAX, true},
        {"BSSN_WAVELET_TOL_FUNCTION_R0", bssn::BSSN_WAVELET_TOL_FUNCTION_R0,
         true},
        {"BSSN_WAVELET_TOL_FUNCTION_R1", bssn::BSSN_WAVELET_TOL_FUNCTION_R1,
         true},
        // these arrays are read in below
        {"BSSN_NUM_REFINE_VARS", bssn::BSSN_NUM_REFINE_VARS, true},
        {"BSSN_NUM_EVOL_VARS_VTU_OUTPUT", bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT,
         true},
        {"BSSN_NUM_CONST_VARS_VTU_OUTPUT", bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT,
         true},
        {"BSSN_CFL_FACTOR", bssn::BSSN_CFL_FACTOR, false},
        {"BSSN_VTU_Z_SLICE_ONLY", bssn::BSSN_VTU_Z_SLICE_ONLY, false},
        // GW_EXTRACT_FREQ is set later because it's optional and needs to be
        // set to something different
        {"BSSN_GW_EXTRACT_FREQ", bssn::BSSN_GW_EXTRACT_FREQ, false},
        {"BSSN_TIME_STEP_OUTPUT_FREQ", bssn::BSSN_TIME_STEP_OUTPUT_FREQ, false},
        {"BSSN_BH1_AMR_R", bssn::BSSN_BH1_AMR_R, false},
        {"BSSN_BH2_AMR_R", bssn::BSSN_BH2_AMR_R, false},
        {"BSSN_AMR_R_RATIO", bssn::BSSN_AMR_R_RATIO, false},
        {"BSSN_DENDRO_AMR_FAC_POST_MERGER",
         bssn::BSSN_DENDRO_AMR_FAC_POST_MERGER, false},
        {"BSSN_BH1_MAX_LEV", bssn::BSSN_BH1_MAX_LEV, false},
        {"BSSN_BH2_MAX_LEV", bssn::BSSN_BH2_MAX_LEV, false},
        {"BSSN_BH2_MAX_LEV", bssn::BSSN_BH2_MAX_LEV, false},
        {"BSSN_INIT_GRID_ITER", bssn::BSSN_INIT_GRID_ITER, false},
        {"BSSN_GW_REFINE_WTOL", bssn::BSSN_GW_REFINE_WTOL, false},
        {"BSSN_MINDEPTH", bssn::BSSN_MINDEPTH, false},
        {"BSSN_BH1_CONSTRAINT_R", bssn::BSSN_BH1_CONSTRAINT_R, false},
        {"BSSN_BH2_CONSTRAINT_R", bssn::BSSN_BH2_CONSTRAINT_R, false},
        {"BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE",
         bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE, false},
        {"BSSN_NYQUIST_M", bssn::BSSN_NYQUIST_M, false},
        {"BSSN_SCALE_VTU_AND_GW_EXTRACTION",
         bssn::BSSN_SCALE_VTU_AND_GW_EXTRACTION, false},
        {"BSSN_SSL_SIGMA", bssn::BSSN_SSL_SIGMA, false},
        {"BSSN_SSL_H", bssn::BSSN_SSL_H, false},
        {"BSSN_EH_COARSEN_VAL", bssn::BSSN_EH_COARSEN_VAL, false},
        {"BSSN_EH_REFINE_VAL", bssn::BSSN_EH_REFINE_VAL, false},

        // TPID parameters
        {"TPID_TARGET_M_PLUS", TPID::target_M_plus, true},
        {"TPID_TARGET_M_MINUS", TPID::target_M_minus, true},
        {"TPID_PAR_B", TPID::par_b, true},
        {"TPID_INITIAL_LAPSE_PSI_EXPONENT", TPID::initial_lapse_psi_exponent,
         true},
        {"TPID_NPOINTS_A", TPID::npoints_A, true},
        {"TPID_NPOINTS_B", TPID::npoints_B, true},
        {"TPID_NPOINTS_PHI", TPID::npoints_phi, true},
        {"TPID_GIVE_BARE_MASS", TPID::give_bare_mass, true},
        {"INITIAL_LAPSE", TPID::initial_lapse, true},
        {"TPID_SOLVE_MOMENTUM_CONSTRAINT", TPID::solve_momentum_constraint,
         true},
        {"TPID_GRID_SETUP_METHOD", TPID::grid_setup_method, true},
        {"TPID_VERBOSE", TPID::verbose, true},
        {"TPID_ADM_TOL", TPID::adm_tol, true},
        {"TPID_NEWTON_TOL", TPID::Newton_tol, true},
        {"TPID_FILEPREFIX", TPID::FILE_PREFIX, false},
        {"TPID_REPLACE_LAPSE_WITH_SQRT_CHI", TPID::replace_lapse_with_sqrt_chi,
         false},

        // wave extraction parameters
        {"EXTRACTION_VAR_ID", BHLOC::EXTRACTION_VAR_ID, false},
        {"EXTRACTION_TOL", BHLOC::EXTRACTION_TOL, false},
        {"BSSN_GW_NUM_RADAII", GW::BSSN_GW_NUM_RADAII, true},
        {"BSSN_GW_NUM_LMODES", GW::BSSN_GW_NUM_LMODES, true},

        // ONE required AH parameter:
        {"AEH_SOLVER_FREQ", AEH::AEH_SOLVER_FREQ, true},
    };

    // then load the parameters
    for (const auto& param : parsList) {
        set_param(parFile, param);
    }

    bssn::BH1 = BH(parFile["BSSN_BH1"]["MASS"].as_floating(),
                   parFile["BSSN_BH1"]["X"].as_floating(),
                   parFile["BSSN_BH1"]["Y"].as_floating(),
                   parFile["BSSN_BH1"]["Z"].as_floating(),
                   parFile["BSSN_BH1"]["V_X"].as_floating(),
                   parFile["BSSN_BH1"]["V_Y"].as_floating(),
                   parFile["BSSN_BH1"]["V_Z"].as_floating(),
                   parFile["BSSN_BH1"]["SPIN"].as_floating(),
                   parFile["BSSN_BH1"]["SPIN_THETA"].as_floating(),
                   parFile["BSSN_BH1"]["SPIN_PHI"].as_floating());
    bssn::BH2 = BH(parFile["BSSN_BH2"]["MASS"].as_floating(),
                   parFile["BSSN_BH2"]["X"].as_floating(),
                   parFile["BSSN_BH2"]["Y"].as_floating(),
                   parFile["BSSN_BH2"]["Z"].as_floating(),
                   parFile["BSSN_BH2"]["V_X"].as_floating(),
                   parFile["BSSN_BH2"]["V_Y"].as_floating(),
                   parFile["BSSN_BH2"]["V_Z"].as_floating(),
                   parFile["BSSN_BH2"]["SPIN"].as_floating(),
                   parFile["BSSN_BH2"]["SPIN_THETA"].as_floating(),
                   parFile["BSSN_BH2"]["SPIN_PHI"].as_floating());
    used_params.insert("BSSN_BH1");
    used_params.insert("BSSN_BH2");

    // update some values based on the ONLY param
    if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY) {
        bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL = false;
    }
    if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL) {
        bssn::BSSN_CAKO_ENABLED = true;
    }

    // force update gw extract frequency if it wasn't actually set!
    if (bssn::BSSN_GW_EXTRACT_FREQ ==
        std::numeric_limits<unsigned int>::max()) {
        bssn::BSSN_GW_EXTRACT_FREQ =
            std::max(1u, bssn::BSSN_IO_OUTPUT_FREQ >> 1u);
    }

    // force update time step output freq if it wasn't set
    if (bssn::BSSN_TIME_STEP_OUTPUT_FREQ ==
        std::numeric_limits<unsigned int>::max()) {
        bssn::BSSN_TIME_STEP_OUTPUT_FREQ = bssn::BSSN_GW_EXTRACT_FREQ;
    }

    // make sure we're limited by the max depth
    if (bssn::BSSN_BH1_MAX_LEV == std::numeric_limits<unsigned int>::max()) {
        bssn::BSSN_BH1_MAX_LEV = bssn::BSSN_MAXDEPTH;
    }
    if (bssn::BSSN_BH2_MAX_LEV == std::numeric_limits<unsigned int>::max()) {
        bssn::BSSN_BH2_MAX_LEV = bssn::BSSN_MAXDEPTH;
    }

    // these ones are a bit trickier than set_param
    bssn::BSSN_LAMBDA[0]   = parFile["BSSN_LAMBDA"][0].as_integer();
    bssn::BSSN_LAMBDA[1]   = parFile["BSSN_LAMBDA"][1].as_integer();
    bssn::BSSN_LAMBDA[2]   = parFile["BSSN_LAMBDA"][2].as_integer();
    bssn::BSSN_LAMBDA[3]   = parFile["BSSN_LAMBDA"][3].as_integer();
    bssn::BSSN_LAMBDA_F[0] = parFile["BSSN_LAMBDA_F"][0].as_floating();
    bssn::BSSN_LAMBDA_F[1] = parFile["BSSN_LAMBDA_F"][1].as_floating();
    used_params.insert("BSSN_LAMBDA");
    used_params.insert("BSSN_LAMBDA_F");

    bssn::BSSN_XI[0] = (unsigned int)parFile["BSSN_XI"][0].as_integer();
    bssn::BSSN_XI[1] = (unsigned int)parFile["BSSN_XI"][1].as_integer();
    bssn::BSSN_XI[2] = (unsigned int)parFile["BSSN_XI"][2].as_integer();
    used_params.insert("BSSN_XI");

    bssn::BSSN_ETA_POWER[0] = parFile["BSSN_ETA_POWER"][0].as_floating();
    bssn::BSSN_ETA_POWER[1] = parFile["BSSN_ETA_POWER"][1].as_floating();
    used_params.insert("BSSN_ETA_POWER");

    for (unsigned int i = 0; i < bssn::BSSN_NUM_REFINE_VARS; i++)
        bssn::BSSN_REFINE_VARIABLE_INDICES[i] =
            parFile["BSSN_REFINE_VARIABLE_INDICES"][i].as_integer();
    used_params.insert("BSSN_REFINE_VARIABLE_INDICES");

    for (unsigned int i = 0; i < bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT; i++)
        bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[i] =
            parFile["BSSN_VTU_OUTPUT_EVOL_INDICES"][i].as_integer();
    used_params.insert("BSSN_VTU_OUTPUT_EVOL_INDICES");

    for (unsigned int i = 0; i < bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT; i++)
        bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i] =
            parFile["BSSN_VTU_OUTPUT_CONST_INDICES"][i].as_integer();
    used_params.insert("BSSN_VTU_OUTPUT_CONST_INDICES");

    bssn::BSSN_IO_OUTPUT_FREQ_TRUE  = bssn::BSSN_IO_OUTPUT_FREQ;
    bssn::BSSN_GW_EXTRACT_FREQ_TRUE = bssn::BSSN_GW_EXTRACT_FREQ;

    /* Setting remainder Parameters for TPID */
    TPID::par_m_plus                = TPID::target_M_plus;
    TPID::par_m_minus               = TPID::target_M_minus;

    TPID::par_P_plus[0] =
        bssn::BH1.getVx();  // parFile["TPID_PAR_P_PLUS"]["X"];
    TPID::par_P_plus[1] =
        bssn::BH1.getVy();  // parFile["TPID_PAR_P_PLUS"]["Y"];
    TPID::par_P_plus[2] =
        bssn::BH1.getVz();  // parFile["TPID_PAR_P_PLUS"]["Z"];

    TPID::par_P_minus[0] =
        bssn::BH2.getVx();  // parFile["TPID_PAR_P_MINUS"]["X"];
    TPID::par_P_minus[1] =
        bssn::BH2.getVy();  // parFile["TPID_PAR_P_MINUS"]["Y"];
    TPID::par_P_minus[2] =
        bssn::BH2.getVz();  // parFile["TPID_PAR_P_MINUS"]["Z"];

    TPID::par_S_plus[0] =
        bssn::BH1.getBHSpin() * sin(bssn::BH1.getBHSpinTheta()) *
        cos(bssn::BH1.getBHSpinPhi());  // parFile["TPID_PAR_S_PLUS"]["X"];
    TPID::par_S_plus[1] =
        bssn::BH1.getBHSpin() * sin(bssn::BH1.getBHSpinTheta()) *
        sin(bssn::BH1.getBHSpinPhi());  // parFile["TPID_PAR_S_PLUS"]["Y"];
    TPID::par_S_plus[2] =
        bssn::BH1.getBHSpin() *
        cos(bssn::BH1.getBHSpinTheta());  // parFile["TPID_PAR_S_PLUS"]["Z"];

    TPID::par_S_minus[0] =
        bssn::BH2.getBHSpin() * sin(bssn::BH2.getBHSpinTheta()) *
        cos(bssn::BH2.getBHSpinPhi());  // parFile["TPID_PAR_S_MINUS"]["X"];
    TPID::par_S_minus[1] =
        bssn::BH2.getBHSpin() * sin(bssn::BH2.getBHSpinTheta()) *
        sin(bssn::BH2.getBHSpinPhi());  // parFile["TPID_PAR_S_MINUS"]["Y"];
    TPID::par_S_minus[2] =
        bssn::BH2.getBHSpin() *
        cos(bssn::BH2.getBHSpinTheta());  // parFile["TPID_PAR_S_MINUS"]["Z"];

    TPID::center_offset[0] = parFile["TPID_CENTER_OFFSET"]["X"].as_floating();
    TPID::center_offset[1] = parFile["TPID_CENTER_OFFSET"]["Y"].as_floating();
    TPID::center_offset[2] = parFile["TPID_CENTER_OFFSET"]["Z"].as_floating();
    used_params.insert("TPID_CENTER_OFFSET");

    // extractio parameters that need to be read separately
    for (unsigned int i = 0; i < GW::BSSN_GW_NUM_RADAII; i++)
        GW::BSSN_GW_RADAII[i] = parFile["BSSN_GW_RADAII"][i].as_floating();
    used_params.insert("BSSN_GW_RADAII");

    for (unsigned int i = 0; i < GW::BSSN_GW_NUM_LMODES; i++)
        GW::BSSN_GW_L_MODES[i] = parFile["BSSN_GW_L_MODES"][i].as_integer();
    used_params.insert("BSSN_GW_L_MODES");

    if (parFile.contains("BSSN_REFINEMENT_MODE"))
        bssn::BSSN_REFINEMENT_MODE = static_cast<bssn::RefinementMode>(
            parFile["BSSN_REFINEMENT_MODE"].as_integer());
    used_params.insert("BSSN_REFINEMENT_MODE");

    BSSN_OCTREE_MAX[0] = (double)(1u << bssn::BSSN_MAXDEPTH);
    BSSN_OCTREE_MAX[1] = (double)(1u << bssn::BSSN_MAXDEPTH);
    BSSN_OCTREE_MAX[2] = (double)(1u << bssn::BSSN_MAXDEPTH);

    BSSN_COMPD_MIN[0]  = bssn::BSSN_GRID_MIN_X;
    BSSN_COMPD_MIN[1]  = bssn::BSSN_GRID_MIN_Y;
    BSSN_COMPD_MIN[2]  = bssn::BSSN_GRID_MIN_Z;

    BSSN_COMPD_MAX[0]  = bssn::BSSN_GRID_MAX_X;
    BSSN_COMPD_MAX[1]  = bssn::BSSN_GRID_MAX_Y;
    BSSN_COMPD_MAX[2]  = bssn::BSSN_GRID_MAX_Z;

    if (BSSN_NUM_REFINE_VARS > BSSN_NUM_VARS) {
        std::cout << "Error[parameter file]: Number of refine variables should "
                     "be less than number of BSSN_NUM_VARS"
                  << std::endl;
        exit(0);
    }
    if (BSSN_NUM_EVOL_VARS_VTU_OUTPUT > BSSN_NUM_VARS) {
        std::cout << "Error[parameter file]: Number of evolution VTU variables "
                     "should be less than number of BSSN_NUM_VARS"
                  << std::endl;
        exit(0);
    }
    if (BSSN_NUM_CONST_VARS_VTU_OUTPUT > BSSN_CONSTRAINT_NUM_VARS) {
        std::cout
            << "Error[parameter file]: Number of constraint VTU variables "
               "should be less than number of BSSN_CONSTRAINT_NUM_VARS"
            << std::endl;
        exit(0);
    }

    BSSN_PADDING_WIDTH = BSSN_ELE_ORDER >> 1u;
    bssn::BSSN_BH_LOC[0] =
        Point(BH1.getBHCoordX(), BH1.getBHCoordY(), BH1.getBHCoordZ());
    bssn::BSSN_BH_LOC[1] =
        Point(BH2.getBHCoordX(), BH2.getBHCoordY(), BH2.getBHCoordZ());
    bssn::BSSN_BH1_MASS = BH1.getBHMass();
    bssn::BSSN_BH2_MASS = BH2.getBHMass();

    // if the parFile has the AEH "dictionary"
    if (parFile.contains("AEH_PARAMS")) {
        auto aeh_pars                                 = parFile["AEH_PARAMS"];

        std::vector<ParameterInformation> aehParsList = {
            {"N_HORIZONS", AEH::N_HORIZONS, false},
            {"N_RESOLUTIONS_MULTIGRID", AEH::N_RESOLUTIONS_MULTIGRID, false},
            {"N_RESOLUTIONS_MULTIGRID", AEH::N_RESOLUTIONS_MULTIGRID, false},

            {"MAX_ITERATIONS", AEH::MAX_ITERATIONS, false},
            {"INITIAL_X_CENTER", AEH::INITIAL_X_CENTER, false},
            {"INITIAL_Y_CENTER", AEH::INITIAL_Y_CENTER, false},
            {"INITIAL_Z_CENTER", AEH::INITIAL_Z_CENTER, false},
            {"M_SCALE", AEH::M_SCALE, false},
            {"CFL_FACTOR", AEH::CFL_FACTOR, false},
            {"THETA_L2_M_TOL", AEH::THETA_L2_M_TOL, false},
            {"THETA_LINF_M_TOL", AEH::THETA_LINF_M_TOL, false},
            {"ETA_DAMP_M", AEH::ETA_DAMP_M, false},
            {"KO_STRENGTH", AEH::KO_STRENGTH, false},
            {"MAX_SEARCH_RADIUS", AEH::MAX_SEARCH_RADIUS, false},
            {"NR_INTERP_MAX", AEH::NR_INTERP_MAX, false},
            {"NPHI_MAX", AEH::NPHI_MAX, false},
            {"NTHETA_MAX", AEH::NTHETA_MAX, false},
            {"AEH_SAVE_DIR", AEH::AEH_SAVE_DIR, false},
            {"NUM_RESOLUTIONS_AFTER_FIND", AEH::NUM_RESOLUTIONS_AFTER_FIND,
             false},
            {"NTHETA_ARRAY", AEH::NTHETA_ARRAY, false},
            {"NPHI_ARRAY", AEH::NPHI_ARRAY, false},
            {"ENABLE_ETA_VARYING_ALG", AEH::ENABLE_ETA_VARYING_ALG, false},
            {"VERBOSITY_LEVEL", AEH::VERBOSITY_LEVEL, false}};

        // then load the parameters
        for (const auto& param : aehParsList) {
            set_param(aeh_pars, param);
        }

        used_params.insert("AEH_PARAMS");
    }

    bool is_bbh_temp = false;
    std::vector<dendro_aeh::SimpleBlackHoleData> simpleBHData;
    if (BSSN_ID_TYPE == 0 || BSSN_ID_TYPE == 1) {
        is_bbh_temp  = true;
        simpleBHData = {dendro_aeh::SimpleBlackHoleData(
                            bssn::BH1.getBHCoordX(), bssn::BH1.getBHCoordY(),
                            bssn::BH1.getBHCoordZ(), bssn::BH1.getBHMass()),
                        dendro_aeh::SimpleBlackHoleData(
                            bssn::BH2.getBHCoordX(), bssn::BH2.getBHCoordY(),
                            bssn::BH2.getBHCoordZ(), bssn::BH2.getBHMass())};
    }

    // this is the function that will convert our inputs of chi, trK, at, and Gt
    // into the standard metric variables of Gd and Kd
    auto transform = [](const std::vector<double>& inputs) {
        // do the transformation here

        // so this is the transformation that'll be applied at each point, it
        // expects 12 outputs in this order: gd00, gd01, gd02, gd11, gd12, gd22,
        // Kd00, Kd01, Kd02, Kd11, Kd12, Kd22
        std::vector<double> output(12, 0.0);

        // the input follow aeh::AEH_INDICES, which are "constant"
        constexpr int chi_idx  = 0;
        constexpr int trK_idx  = 1;
        constexpr int at00_idx = 2;
        constexpr int at01_idx = 3;
        constexpr int at02_idx = 4;
        constexpr int at11_idx = 5;
        constexpr int at12_idx = 6;
        constexpr int at22_idx = 7;
        constexpr int Gt00_idx = 8;
        constexpr int Gt01_idx = 9;
        constexpr int Gt02_idx = 10;
        constexpr int Gt11_idx = 11;
        constexpr int Gt12_idx = 12;
        constexpr int Gt22_idx = 13;

        constexpr int gd00     = 0;
        constexpr int gd01     = 1;
        constexpr int gd02     = 2;
        constexpr int gd11     = 3;
        constexpr int gd12     = 4;
        constexpr int gd22     = 5;
        constexpr int Kd00     = 6;
        constexpr int Kd01     = 7;
        constexpr int Kd02     = 8;
        constexpr int Kd11     = 9;
        constexpr int Kd12     = 10;
        constexpr int Kd22     = 11;

        // chi = idetgd^(1/3)
        double idetgd          = pow(inputs[chi_idx], -3.0);
        // detgd = 1 / idetgd
        double detgd           = 1.0 / idetgd;

        double inv_chi         = 1.0 / inputs[chi_idx];
        double trK_third       = (1.0 / 3.0) * inputs[trK_idx];

        // calculate local gd
        output[gd00]           = inputs[Gt00_idx] * inv_chi;
        output[gd01]           = inputs[Gt01_idx] * inv_chi;
        output[gd02]           = inputs[Gt02_idx] * inv_chi;
        output[gd11]           = inputs[Gt11_idx] * inv_chi;
        output[gd12]           = inputs[Gt12_idx] * inv_chi;
        output[gd22]           = inputs[Gt22_idx] * inv_chi;

        // fill in Kd values
        // basically: Atd[i,j] / chi + (1/3) * gd[i,j] * trK
        output[Kd00] = inputs[at00_idx] * inv_chi + output[gd00] * trK_third;
        output[Kd01] = inputs[at01_idx] * inv_chi + output[gd01] * trK_third;
        output[Kd02] = inputs[at02_idx] * inv_chi + output[gd02] * trK_third;
        output[Kd11] = inputs[at11_idx] * inv_chi + output[gd11] * trK_third;
        output[Kd12] = inputs[at12_idx] * inv_chi + output[gd12] * trK_third;
        output[Kd22] = inputs[at22_idx] * inv_chi + output[gd22] * trK_third;

        return output;
    };

    Point grid_limits[2];
    Point domain_limits[2];

    grid_limits[0]   = Point(bssn::BSSN_OCTREE_MIN[0], bssn::BSSN_OCTREE_MIN[1],
                             bssn::BSSN_OCTREE_MIN[2]);
    grid_limits[1]   = Point(bssn::BSSN_OCTREE_MAX[0], bssn::BSSN_OCTREE_MAX[1],
                             bssn::BSSN_OCTREE_MAX[2]);

    domain_limits[0] = Point(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                             bssn::BSSN_COMPD_MIN[2]);
    domain_limits[1] = Point(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                             bssn::BSSN_COMPD_MAX[2]);

    // now we build up the aeh solver
    AEH::aeh         = std::make_unique<dendro_aeh::AEH_BHaHAHA>(
        AEH::N_HORIZONS, is_bbh_temp, AEH::INITIAL_X_CENTER,
        AEH::INITIAL_Y_CENTER, AEH::INITIAL_Z_CENTER,
        AEH::N_RESOLUTIONS_MULTIGRID, AEH::M_SCALE, AEH::CFL_FACTOR,
        AEH::MAX_ITERATIONS, AEH::THETA_L2_M_TOL, AEH::THETA_LINF_M_TOL,
        AEH::ETA_DAMP_M, AEH::KO_STRENGTH, AEH::MAX_SEARCH_RADIUS,
        AEH::NR_INTERP_MAX, AEH::NTHETA_MAX, AEH::NPHI_MAX, AEH::AEH_SAVE_DIR,
        simpleBHData, AEH::AEH_INDICES, transform, grid_limits, domain_limits,
        AEH::NUM_RESOLUTIONS_AFTER_FIND, AEH::NTHETA_ARRAY, AEH::NPHI_ARRAY,
        AEH::ENABLE_ETA_VARYING_ALG, AEH::VERBOSITY_LEVEL);

    // check on the unused parameters to see what we didn't use
    std::unordered_set<std::string> all_params_in_parfile;
    for (const auto& [key, _] : parFile.as_table()) {
        all_params_in_parfile.insert(key);
    }

    // find the unused ones
    std::vector<std::string> unused_params;
    for (const auto& key : all_params_in_parfile) {
        if (used_params.find(key) == used_params.end()) {
            unused_params.push_back(key);
        }
    }

    // warn about parameters that are in the TOML file that we don't actually
    // need
    if (!unused_params.empty() && rank == 0) {
        std::cout << YLW
                  << "WARNING: Unused parameters in TOML file:" << std::endl;

        for (const auto& key : unused_params) {
            std::cout << "\t -- " << key << std::endl;
        }
        std::cout << NRM << std::endl;
    }

    MPI_Barrier(comm);
}

}  // namespace bssn

namespace TPID {
double target_M_plus              = 1.0;
double target_M_minus             = 1.0;
double par_m_plus                 = 1.0;
double par_m_minus                = 1.0;
double par_b                      = 4.0;
double par_P_plus[3]              = {0.0, 0.0, 0.0};
double par_P_minus[3]             = {0.0, 0.0, 0.0};
double par_S_plus[3]              = {0.0, 0.0, 0.0};
double par_S_minus[3]             = {0.0, 0.0, 0.0};
double center_offset[3]           = {0.0, 0.0, 0.00014142135623730951};
double initial_lapse_psi_exponent = -2;
int npoints_A                     = 30;
int npoints_B                     = 30;
int npoints_phi                   = 16;
int give_bare_mass                = 0;
int initial_lapse                 = 2;
int solve_momentum_constraint     = 1;
int grid_setup_method             = 1;
int verbose                       = 1;
double adm_tol                    = 1.0e-10;
double Newton_tol                 = 1.0e-10;
std::string FILE_PREFIX           = "tpid";
bool replace_lapse_with_sqrt_chi  = false;
}  // namespace TPID

namespace BHLOC {
unsigned int EXTRACTION_VAR_ID = bssn::VAR::U_ALPHA;
double EXTRACTION_TOL          = 0.3;
}  // namespace BHLOC

namespace GW {

unsigned int BSSN_GW_NUM_RADAII;

unsigned int BSSN_GW_NUM_LMODES;

double BSSN_GW_RADAII[BSSN_GW_MAX_RADAII];

unsigned int BSSN_GW_SPIN = 2;

unsigned int BSSN_GW_L_MODES[BSSN_GW_MAX_LMODES];
}  // namespace GW

namespace AEH {

std::unique_ptr<dendro_aeh::AEH_BHaHAHA> aeh = nullptr;

unsigned int N_HORIZONS                      = 3;
unsigned int N_RESOLUTIONS_MULTIGRID         = 3;
std::vector<int> MAX_ITERATIONS              = {10000, 10000, 10000};
std::vector<double> INITIAL_X_CENTER         = {0.0, 0.0, 0.0};
std::vector<double> INITIAL_Y_CENTER         = {0.0, 0.0, 0.0};
std::vector<double> INITIAL_Z_CENTER         = {0.0, 0.0, 0.0};
std::vector<double> M_SCALE                  = {0.0, 0.0, 0.0};
std::vector<double> CFL_FACTOR               = {1.0, 1.0, 1.0};
std::vector<double> THETA_L2_M_TOL           = {1.e-5, 1.e-5, 1.e-5};
std::vector<double> THETA_LINF_M_TOL         = {1.e-2, 1.e-2, 1.e-2};
std::vector<double> ETA_DAMP_M               = {7.0, 7.0, 7.0};
std::vector<double> KO_STRENGTH              = {0.0, 0.0, 0.0};
std::vector<double> MAX_SEARCH_RADIUS        = {1.5, 1.5, 1.5};
std::vector<int> NR_INTERP_MAX               = {48, 48, 48};
int NTHETA_MAX                               = 32;
int NPHI_MAX                                 = 64;
std::string AEH_SAVE_DIR                     = "";
int NUM_RESOLUTIONS_AFTER_FIND               = 3;
std::vector<int> NTHETA_ARRAY                = {8, 16, 32};
std::vector<int> NPHI_ARRAY                  = {16, 32, 64};
int ENABLE_ETA_VARYING_ALG                   = 0;
int VERBOSITY_LEVEL                          = 1;

unsigned int AEH_LMAX                        = 6;
unsigned int AEH_Q_THETA                     = 32;
unsigned int AEH_Q_PHI                       = 32;
unsigned int AEH_MAXITER                     = 50;
double AEH_ATOL                              = 1e-8;
double AEH_RTOL                              = 1e-8;
unsigned int AEH_SOLVER_FREQ                 = 0;

double AEH_ALPHA                             = 1.0;
double AEH_BETA                              = 0.1;

}  // namespace AEH
