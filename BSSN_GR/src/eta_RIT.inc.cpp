// RIT constant eta formulation
const double w        = r_coord / bssn::RIT_ETA_WIDTH;
const double arg      = -w * w * w * w;
const double eta_rit  = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER) * exp(arg) + bssn::RIT_ETA_OUTER;
const double eta      = (bssn::BSSN_CCZ4_ETA >= 0.0) ? bssn::BSSN_CCZ4_ETA : eta_rit;
